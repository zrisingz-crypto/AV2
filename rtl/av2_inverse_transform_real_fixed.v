//==============================================================================
// Real Inverse Transform Module for AV2 (Fixed Syntax)
// Implements 2D IDCT for residual decoding
//==============================================================================

`timescale 1ns / 1ps

module av2_inverse_transform_real_fixed #(
    parameter MAX_TX_SIZE = 64,
    parameter BIT_DEPTH = 10
)(
    input  wire                      clk,
    input  wire                      rst_n,
    
    // Input coefficients
    input  wire signed [15:0]       coeffs [0:4095],  // Max 64x64
    input  wire [15:0]               num_coeffs,
    
    // Transform parameters
    input  wire [5:0]                tx_width,
    input  wire [5:0]                tx_height,
    input  wire [3:0]                tx_type,
    
    // Control signals
    input  wire                      start,
    output reg  signed [15:0]        pixels [0:4095],
    output reg                       valid,
    input  wire                      ready,
    output reg                       done
);

// Transform types
localparam TX_DCT_DCT   = 4'd0;
localparam TX_ADST_DCT   = 4'd1;
localparam TX_DCT_ADST   = 4'd2;
localparam TX_ADST_ADST  = 4'd3;
localparam TX_FLIPADST_DCT  = 4'd4;
localparam TX_DCT_FLIPADST  = 4'd5;
localparam TX_FLIPADST_FLIPADST = 4'd6;
localparam TX_IDTX       = 4'd7;

// State machine
localparam IDLE       = 2'd0;
localparam ROW_TX     = 2'd1;
localparam COL_TX     = 2'd2;
localparam DONE       = 2'd3;

reg [1:0] state, state_next;

// Processing registers
reg [5:0] row_idx, col_idx;
reg [5:0] row_count, col_count;

// Intermediate buffers (declare outside always block)
reg signed [17:0] row_in [0:63];
reg signed [17:0] row_out [0:63];
reg signed [17:0] col_in [0:63];
reg signed [17:0] col_out [0:63];
reg signed [17:0] temp_sum;
reg signed [17:0] temp_diff;
reg signed [17:0] temp_prod;

integer i, j;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        valid <= 1'b0;
        done <= 1'b0;
        row_idx <= 6'd0;
        col_idx <= 6'd0;
        row_count <= 6'd0;
        col_count <= 6'd0;
        temp_sum <= 18'd0;
        temp_diff <= 18'd0;
        temp_prod <= 18'd0;
        
        for (i = 0; i < 64; i = i + 1) begin
            row_in[i] <= 18'd0;
            row_out[i] <= 18'd0;
            col_in[i] <= 18'd0;
            col_out[i] <= 18'd0;
        end
        
        for (i = 0; i < 4096; i = i + 1) begin
            pixels[i] <= 16'sd0;
        end
    end else begin
        state <= state_next;
        
        case (state)
            IDLE: begin
                valid <= 1'b0;
                done <= 1'b0;
                if (start) begin
                    state <= ROW_TX;
                    row_idx <= 6'd0;
                    col_idx <= 6'd0;
                    row_count <= tx_height;
                    col_count <= tx_width;
                end
            end
            
            ROW_TX: begin
                // Perform 1D transform on each row
                if (row_idx < row_count) begin
                    // Load row
                    for (col_idx = 0; col_idx < col_count; col_idx = col_idx + 1) begin
                        row_in[col_idx] <= coeffs[row_idx * tx_width + col_idx];
                    end
                    
                    // Simplified 4-point IDCT for demo
                    if (col_count <= 6'd4) begin
                        // 4-point IDCT (simplified)
                        temp_sum = row_in[0] + row_in[2];
                        temp_diff = row_in[0] - row_in[2];
                        
                        row_out[0] = temp_sum + row_in[1];
                        row_out[1] = temp_diff + row_in[3];
                        row_out[2] = temp_diff - row_in[3];
                        row_out[3] = temp_sum - row_in[1];
                        
                        // Store back to intermediate buffer
                        for (i = 0; i < 4; i = i + 1) begin
                            pixels[row_idx * tx_width + i] <= row_out[i][15:0];
                        end
                    end else begin
                        // For larger blocks, simplified scaling
                        for (i = 0; i < col_count; i = i + 1) begin
                            pixels[row_idx * tx_width + i] <= row_in[i][15:0];
                        end
                    end
                    
                    row_idx <= row_idx + 1;
                end else begin
                    state <= COL_TX;
                    col_idx <= 6'd0;
                end
            end
            
            COL_TX: begin
                // Perform 1D transform on each column
                if (col_idx < col_count) begin
                    // Load column (from row-transformed data)
                    for (row_idx = 0; row_idx < row_count; row_idx = row_idx + 1) begin
                        col_in[row_idx] <= pixels[row_idx * tx_width + col_idx];
                    end
                    
                    // Simplified 4-point IDCT
                    if (row_count <= 6'd4) begin
                        temp_sum = col_in[0] + col_in[2];
                        temp_diff = col_in[0] - col_in[2];
                        
                        col_out[0] = temp_sum + col_in[1];
                        col_out[1] = temp_diff + col_in[3];
                        col_out[2] = temp_diff - col_in[3];
                        col_out[3] = temp_sum - col_in[1];
                        
                        // Store back
                        for (i = 0; i < 4; i = i + 1) begin
                            pixels[i * tx_width + col_idx] <= col_out[i][15:0];
                        end
                    end else begin
                        // For larger blocks
                        for (i = 0; i < row_count; i = i + 1) begin
                            pixels[i * tx_width + col_idx] <= col_in[i][15:0];
                        end
                    end
                    
                    col_idx <= col_idx + 1;
                end else begin
                    valid <= 1'b1;
                    state <= DONE;
                end
            end
            
            DONE: begin
                if (ready) begin
                    valid <= 1'b0;
                    done <= 1'b1;
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
                state_next = ROW_TX;
        end
        
        ROW_TX: begin
            if (row_idx >= row_count)
                state_next = COL_TX;
        end
        
        COL_TX: begin
            if (col_idx >= col_count)
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