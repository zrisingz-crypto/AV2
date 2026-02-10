//==============================================================================
// Inverse Transform Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_inverse_transform #(
    parameter MAX_TX_SIZE = 64
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire signed [15:0]      coeff_in[0:4095],
    input  wire [5:0]             tx_width,
    input  wire [5:0]             tx_height,
    input  wire [3:0]             tx_type,
    input  wire                    start,
    output reg  signed [15:0]      pixel_out[0:4095],
    output reg                     valid,
    input  wire                    ready
);

// Simple stub: pass-through coefficients as pixels
reg busy;
reg [3:0] cycle_count;

integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        valid <= 1'b0;
        busy <= 1'b0;
        cycle_count <= 4'd0;
        for (i = 0; i < 4096; i = i + 1) begin
            pixel_out[i] <= 16'sd0;
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