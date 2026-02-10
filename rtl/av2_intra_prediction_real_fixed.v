//==============================================================================
// Real Intra Prediction Module for AV2 (Fixed Syntax)
// Supports: DC, V, H, PAETH, SMOOTH, SMOOTH_V, SMOOTH_H, ANGULAR
//==============================================================================

`timescale 1ns / 1ps

module av2_intra_prediction_real_fixed #(
    parameter MAX_BLOCK_SIZE = 64,
    parameter BIT_DEPTH = 10
)(
    input  wire                      clk,
    input  wire                      rst_n,
    
    // Reference pixels
    input  wire [9:0]               ref_top [0:MAX_BLOCK_SIZE-1],
    input  wire [9:0]               ref_left [0:MAX_BLOCK_SIZE-1],
    input  wire [9:0]               ref_top_left,
    
    // Prediction mode
    input  wire [6:0]                intra_mode,
    
    // Block dimensions
    input  wire [5:0]                block_width,
    input  wire [5:0]                block_height,
    
    // Control signals
    input  wire                      start,
    output reg  [9:0]               pred_pixels [0:4095],
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
reg [19:0] dc_sum;
reg [9:0]  dc_value;
reg [5:0]  pixel_count;

// Intermediate calculation registers (declare outside always block)
reg [10:0] temp_pred_11;
reg [19:0] smooth_sum;
reg [11:0] smooth_pred;
reg [10:0] base_value;
reg [10:0] p_left, p_top, p_topleft;
reg [6:0]  w_h, w_v;
reg [10:0] pred_h, pred_v;
reg [9:0]  below_pred;
reg [9:0]  right_pred;
reg [5:0]  idx;
reg [10:0] temp_recon;

integer i, j;

reg [31:0] debug_cycle;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n)
        debug_cycle <= 0;
    else
        debug_cycle <= debug_cycle + 1;
end

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        valid <= 1'b0;
        row <= 6'd0;
        col <= 6'd0;
        dc_sum <= 20'd0;
        dc_value <= 10'd0;
        pixel_count <= 6'd0;
        temp_pred_11 <= 11'd0;
        smooth_sum <= 20'd0;
        smooth_pred <= 12'd0;
        base_value <= 11'd0;
        p_left <= 11'd0;
        p_top <= 11'd0;
        p_topleft <= 11'd0;
        w_h <= 7'd0;
        w_v <= 7'd0;
        pred_h <= 11'd0;
        pred_v <= 11'd0;
        below_pred <= 10'd0;
        right_pred <= 10'd0;
        idx <= 6'd0;
        temp_recon <= 11'd0;
        
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
                    if (debug_cycle < 200)
                        $display("[TIME %0t] INTRA_PRED: IDLE -> PREDICTING, ready=%0b", $time, ready);
                end
            end
            
            PREDICTING: begin
                if (debug_cycle < 200 && row == 0 && col == 0)
                    $display("[TIME %0t] INTRA_PRED: PREDICTING, ready=%0b, mode=%0d", $time, ready, intra_mode);
                
                if (ready) begin
                    case (intra_mode)
                        MODE_DC: begin
                            if (row == 6'd0 && col == 6'd0) begin
                                // Calculate DC value
                                dc_sum <= 20'd0;
                                for (i = 0; i < block_width; i = i + 1) begin
                                    dc_sum <= dc_sum + ref_top[i];
                                end
                                for (j = 0; j < block_height; j = j + 1) begin
                                    dc_sum <= dc_sum + ref_left[j];
                                end
                                dc_value <= dc_sum / (block_width + block_height);
                                // Write first pixel immediately with default 128
                                // (will be overwritten in next cycle with correct dc_value)
                                pred_pixels[0] <= 10'd128;
                                col <= col + 1;
                            end else if (row < block_height) begin
                                pred_pixels[row * block_width + col] <= dc_value;
                                
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
                        
                        MODE_V: begin
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
                            if (row < block_height) begin
                                temp_pred_11 = ref_top[col] + ref_left[row] - ref_top_left;
                                
                                p_left   = (temp_pred_11 >= ref_left[row]) ? 
                                           (temp_pred_11 - ref_left[row]) : 
                                           (ref_left[row] - temp_pred_11);
                                p_top     = (temp_pred_11 >= ref_top[col]) ? 
                                           (temp_pred_11 - ref_top[col]) : 
                                           (ref_top[col] - temp_pred_11);
                                p_topleft = (temp_pred_11 >= ref_top_left) ? 
                                           (temp_pred_11 - ref_top_left) : 
                                           (ref_top_left - temp_pred_11);
                                
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
                            if (row < block_height) begin
                                w_h = (col * 32) / block_width;
                                w_v = (row * 32) / block_height;
                                
                                pred_h = ref_left[row] + ((ref_top_left - ref_left[row]) * w_h) >> 5;
                                pred_v = ref_top[col] + ((ref_top_left - ref_top[col]) * w_v) >> 5;
                                
                                smooth_sum = pred_h + pred_v;
                                smooth_pred = smooth_sum >> 1;
                                
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
                            if (row < block_height) begin
                                w_v = (row * 32) / block_height;
                                below_pred = ref_left[block_height - 1];
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
                            if (row < block_height) begin
                                w_h = (col * 32) / block_width;
                                right_pred = ref_top[block_width - 1];
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
                            if (row < block_height) begin
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
                if (debug_cycle < 200)
                    $display("[TIME %0t] INTRA_PRED: DONE, ready=%0b", $time, ready);
                if (ready) begin
                    valid <= 1'b0;
                    state <= IDLE;
                end
            end
        endcase
    end
end

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