//==============================================================================
// CDEF Filter Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_cdef_filter #(
    parameter BLOCK_SIZE = 8
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [9:0]             src_block[0:63],
    input  wire [2:0]             strength_y,
    input  wire [2:0]             strength_uv,
    input  wire [2:0]             damping,
    input  wire                    is_chroma,
    input  wire                    start,
    output reg  [9:0]             dst_block[0:63],
    output reg                     valid,
    input  wire                    ready
);

// Simple stub: pass-through
reg busy;
reg [3:0] cycle_count;

integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        valid <= 1'b0;
        busy <= 1'b0;
        cycle_count <= 4'd0;
        for (i = 0; i < 64; i = i + 1) begin
            dst_block[i] <= 10'd0;
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