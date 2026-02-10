//==============================================================================
// Real Intra Prediction Module for AV2
// Supports: DC, V, H, PAETH, SMOOTH, SMOOTH_V, SMOOTH_H, ANGULAR
//==============================================================================

`timescale 1ns / 1ps

module av2_intra_prediction_real #(
    parameter MAX_BLOCK_SIZE = 64,
    parameter BIT_DEPTH = 10
)(
    input  wire                      clk,
    input  wire                      rst_n,
    
    // Reference pixels
    input  wire [9:0]               ref_top [0:MAX_BLOCK_SIZE-1],      // Above pixels
    input  wire [9:0]               ref_left [0:MAX_BLOCK_SIZE-1],     // Left pixels
    input  wire [9:0]               ref_top_left,                          // Top-left pixel
    
    // Prediction mode
    input  wire [6:0]                intra_mode,      // 0:DC, 1:V, 2:H, 3:PAETH, 
                                                  // 4:SMOOTH, 5:SMOOTH_V, 6:SMOOTH_H,
                                                  // 7-10:ANGULAR (DIAG_45, etc.)
    
    // Block dimensions
    input  wire [5:0]                block_width,
    input  wire [5:0]                block_height,
    
    // Control signals
    input  wire                      start,
    output reg  [9:0]               pred_pixels [0:4095],  // Max 64x64
    output reg                       valid,
    input  wire                      ready
);

// Prediction mode encoding
localparam MODE_DC          = 7'd0;
localparam MODE_V           = 7'd1;
localparam MODE_H           = 7'd2;
localparam MODE_PAETH       = 7'd3;
localparam MODE_SMOOTH      = 7'd4;
localparam MODE_SMOOTH_V    = 7'd5;
localparam MODE_SMOOTH_H    = 7'd6;
localparam MODE_ANGULAR     = 7'd7;

// State machine
localparam IDLE       = 2'd0;
localparam PREDICTING = 2'd1;
localparam DONE       = 2'd2;

reg [1:0] state, state_next;

// Processing registers
reg [5:0] row, col;
reg [19:0] dc_sum;       // Sum for DC prediction (20 bits for safety)
reg [9:0]  dc_value;     // Final DC value
reg [5:0]  pixel_count;  // Counter for pixels processed

// Intermediate calculation registers
reg [10:0] temp_pred_11;  // 11-bit for intermediate calculations
reg [19:0] smooth_sum;    // 20-bit for smooth prediction
reg [11:0] smooth_pred;    // 12-bit intermediate

// Paeth prediction
reg [10:0] base_value;
reg [10:0] p_left, p_top, p_topleft;

integer i, j;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        valid <= 1'b0;
        row <= 6'd0;
        col <= 6'd0;
        dc_sum <= 20'd0;
        dc_value <= 10'd0;
        pixel_count <= 6'd0;
        for (i = 0; i < 4096; i = i + 1) begin
            pred_pixels[i] <= 10'd0;
        end
    end else begin
        state <= state_next;
        
        case (state)
            IDLE: begin
                valid <= 1'b0;
                if (start) begin
                    state <= PREDICTING;
                    row <= 6'd0;
                    col <= 6'd0;
                    dc_sum <= 20'd0;
                    pixel_count <= 6'd0;
                end
            end
            
            PREDICTING: begin
                if (ready) begin
                    // Process prediction based on mode
                    case (intra_mode)
                        MODE_DC: begin
                            // DC prediction: average of top and left reference pixels
                            if (row == 6'd0 && col == 6'd0) begin
                                // First cycle: calculate DC value
                                // Sum top reference pixels
                                for (i = 0; i < block_width; i = i + 1) begin
                                    dc_sum <= dc_sum + ref_top[i];
                                end
                                // Sum left reference pixels
                                for (j = 0; j < block_height; j = j + 1) begin
                                    dc_sum <= dc_sum + ref_left[j];
                                end
                                // Calculate average with rounding
                                dc_value <= dc_sum[19:0] / (block_width + block_height);
                            end else if (row < block_height) begin
                                // Fill block with DC value
                                pred_pixels[row * block_width + col] <= dc_value;
                                
                                // Increment counters
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                // Done
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        MODE_V: begin
                            // Vertical prediction: copy from top reference
                            if (row < block_height) begin
                                pred_pixels[row * block_width + col] <= ref_top[col];
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        MODE_H: begin
                            // Horizontal prediction: copy from left reference
                            if (row < block_height) begin
                                pred_pixels[row * block_width + col] <= ref_left[row];
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        MODE_PAETH: begin
                            // Paeth prediction: base = top + left - topleft
                            // Choose nearest to base among top, left, topleft
                            if (row < block_height) begin
                                temp_pred_11 = ref_top[col] + ref_left[row] - ref_top_left;
                                
                                // Calculate absolute differences
                                p_left   = (temp_pred_11 >= ref_left[row]) ? 
                                           (temp_pred_11 - ref_left[row]) : 
                                           (ref_left[row] - temp_pred_11);
                                p_top     = (temp_pred_11 >= ref_top[col]) ? 
                                           (temp_pred_11 - ref_top[col]) : 
                                           (ref_top[col] - temp_pred_11);
                                p_topleft = (temp_pred_11 >= ref_top_left) ? 
                                           (temp_pred_11 - ref_top_left) : 
                                           (ref_top_left - temp_pred_11);
                                
                                // Choose nearest
                                if (p_left <= p_top && p_left <= p_topleft) begin
                                    pred_pixels[row * block_width + col] <= ref_left[row];
                                end else if (p_top <= p_topleft) begin
                                    pred_pixels[row * block_width + col] <= ref_top[col];
                                end else begin
                                    pred_pixels[row * block_width + col] <= ref_top_left;
                                end
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        MODE_SMOOTH: begin
                            // Smooth prediction: blend from 4 corners
                            // Simplified version for RTL implementation
                            if (row < block_height) begin
                                // Simple bilinear interpolation
                                // Horizontal weight based on column position
                                // Vertical weight based on row position
                                reg [6:0] w_h, w_v;
                                w_h = (col * 32) / block_width;      // 0-32
                                w_v = (row * 32) / block_height;     // 0-32
                                
                                reg [10:0] pred_h, pred_v;
                                pred_h = ref_left[row] + ((ref_top_left - ref_left[row]) * w_h) >> 5;
                                pred_v = ref_top[col] + ((ref_top_left - ref_top[col]) * w_v) >> 5;
                                
                                smooth_sum = pred_h + pred_v;
                                smooth_pred = smooth_sum[19:0] >> 1;  // Average
                                
                                // Clamp to bit depth
                                if (smooth_pred > (10'd1 << BIT_DEPTH) - 1) begin
                                    pred_pixels[row * block_width + col] <= (10'd1 << BIT_DEPTH) - 1;
                                end else begin
                                    pred_pixels[row * block_width + col] <= smooth_pred[9:0];
                                end
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        MODE_SMOOTH_V: begin
                            // Smooth vertical: blend top and bottom-left
                            if (row < block_height) begin
                                reg [6:0] w_v;
                                w_v = (row * 32) / block_height;  // 0-32
                                
                                reg [9:0] below_pred;
                                below_pred = ref_left[block_height - 1];  // Bottom-left reference
                                
                                reg [10:0] pred_v;
                                pred_v = below_pred + ((ref_top[col] - below_pred) * w_v) >> 5;
                                
                                if (pred_v > (10'd1 << BIT_DEPTH) - 1) begin
                                    pred_pixels[row * block_width + col] <= (10'd1 << BIT_DEPTH) - 1;
                                end else begin
                                    pred_pixels[row * block_width + col] <= pred_v[9:0];
                                end
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        MODE_SMOOTH_H: begin
                            // Smooth horizontal: blend left and top-right
                            if (row < block_height) begin
                                reg [6:0] w_h;
                                w_h = (col * 32) / block_width;  // 0-32
                                
                                reg [9:0] right_pred;
                                right_pred = ref_top[block_width - 1];  // Top-right reference
                                
                                reg [10:0] pred_h;
                                pred_h = right_pred + ((ref_left[row] - right_pred) * w_h) >> 5;
                                
                                if (pred_h > (10'd1 << BIT_DEPTH) - 1) begin
                                    pred_pixels[row * block_width + col] <= (10'd1 << BIT_DEPTH) - 1;
                                end else begin
                                    pred_pixels[row * block_width + col] <= pred_h[9:0];
                                end
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        MODE_ANGULAR: begin
                            // Angular prediction (simplified: diagonal 45)
                            // For full implementation, need to support all angles
                            if (row < block_height) begin
                                // Diagonal 45: pred[r][c] = ref_top[c+1] if r < c else ref_left[r-c-1]
                                reg [5:0] idx;
                                if (row <= col) begin
                                    idx = col - row;
                                    if (idx < block_width) begin
                                        pred_pixels[row * block_width + col] <= ref_top[idx];
                                    end else begin
                                        pred_pixels[row * block_width + col] <= ref_top_left;
                                    end
                                end else begin
                                    idx = row - col - 1;
                                    if (idx < block_height) begin
                                        pred_pixels[row * block_width + col] <= ref_left[idx];
                                    end else begin
                                        pred_pixels[row * block_width + col] <= ref_left[block_height - 1];
                                    end
                                end
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                        
                        default: begin
                            // Default to DC prediction
                            if (row < block_height) begin
                                pred_pixels[row * block_width + col] <= 10'd128;
                                
                                if (col < block_width - 1) begin
                                    col <= col + 1;
                                end else begin
                                    col <= 6'd0;
                                    row <= row + 1;
                                end
                            end else begin
                                valid <= 1'b1;
                                state <= DONE;
                            end
                        end
                    endcase
                end
            end
            
            DONE: begin
                if (ready) begin
                    valid <= 1'b0;
                    state <= IDLE;
                end
            end
        endcase
    end
end

// State transition logic
always @(*) begin
    state_next = state;
    
    case (state)
        IDLE: begin
            if (start)
                state_next = PREDICTING;
        end
        
        PREDICTING: begin
            if (valid)
                state_next = DONE;
        end
        
        DONE: begin
            if (ready)
                state_next = IDLE;
        end
        
        default: begin
            state_next = IDLE;
        end
    endcase
end

endmodule