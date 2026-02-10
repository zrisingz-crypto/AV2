//==============================================================================
// Simplified Inverse Transform Module for testing
//==============================================================================

`timescale 1ns / 1ps

module av2_inverse_transform_simple #(
    parameter MAX_TX_SIZE = 64,
    parameter BIT_DEPTH = 10
)(
    input  wire                      clk,
    input  wire                      rst_n,
    
    // Input coefficients
    input  wire signed [15:0]       coeffs [0:4095],
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

// State machine
localparam IDLE     = 2'd0;
localparam PROCESS  = 2'd1;
localparam DONE_S   = 2'd2;

reg [1:0] state;
reg [5:0] row_idx;
reg [5:0] row_count;

integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        valid <= 1'b0;
        done <= 1'b0;
        row_idx <= 6'd0;
        
        for (i = 0; i < 4096; i = i + 1)
            pixels[i] <= 16'sd0;
    end else begin
        case (state)
            IDLE: begin
                valid <= 1'b0;
                done <= 1'b0;
                row_idx <= 6'd0;
                
                if (start) begin
                    row_count <= tx_height;
                    state <= PROCESS;
                end
            end
            
            PROCESS: begin
                // Simple pass-through: copy coeffs to pixels with scaling
                if (row_idx < row_count) begin
                    for (i = 0; i < tx_width; i = i + 1) begin
                        // Scale coefficient to pixel value (simulated)
                        if (coeffs[row_idx * tx_width + i] != 16'sdx)
                            pixels[row_idx * tx_width + i] <= coeffs[row_idx * tx_width + i] * 10;
                        else
                            pixels[row_idx * tx_width + i] <= 16'sd10;  // Default value
                    end
                    row_idx <= row_idx + 1;
                end else begin
                    valid <= 1'b1;
                    state <= DONE_S;
                end
            end
            
            DONE_S: begin
                if (ready) begin
                    valid <= 1'b0;
                    done <= 1'b1;
                    state <= IDLE;
                end
            end
        endcase
    end
end

endmodule
