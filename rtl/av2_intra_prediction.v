//==============================================================================
// Intra Prediction Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_intra_prediction #(
    parameter MAX_BLOCK_SIZE = 64
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [9:0]             ref_top[0:127],
    input  wire [9:0]             ref_left[0:127],
    input  wire [9:0]             ref_top_left,
    input  wire [6:0]             intra_mode,
    input  wire [5:0]             block_width,
    input  wire [5:0]             block_height,
    input  wire                    start,
    output reg  [9:0]             pred_pixels[0:4095],
    output reg                     valid,
    input  wire                    ready
);

// Simple stub: DC prediction (all pixels = 128)
reg busy;
reg [3:0] cycle_count;

integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        valid <= 1'b0;
        busy <= 1'b0;
        cycle_count <= 4'd0;
        for (i = 0; i < 4096; i = i + 1) begin
            pred_pixels[i] <= 10'd128;
        end
    end else begin
        if (start && !busy) begin
            busy <= 1'b1;
            cycle_count <= 4'd0;
        end else if (busy) begin
            if (cycle_count < 4'd4) begin
                cycle_count <= cycle_count + 1;
            end else begin
                valid <= 1'b1;
                if (ready) begin
                    valid <= 1'b0;
                    busy <= 1'b0;
                end
            end
        end
    end
end

endmodule