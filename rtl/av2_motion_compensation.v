//==============================================================================
// Motion Compensation Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_motion_compensation #(
    parameter MAX_WIDTH = 128,
    parameter MAX_HEIGHT = 128,
    parameter MAX_BLOCK_SIZE = 128
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [31:0]             ref_frame_addr,
    output reg  [31:0]             ref_read_addr,
    input  wire [9:0]             ref_pixel_data,
    output reg                      ref_read_en,
    input  wire signed [15:0]      mv_x,
    input  wire signed [15:0]      mv_y,
    input  wire [6:0]              block_width,
    input  wire [6:0]              block_height,
    input  wire [15:0]             block_x,
    input  wire [15:0]             block_y,
    input  wire [3:0]              interp_filter,
    input  wire                    start,
    input  wire                    use_bidir,
    output reg  [9:0]             pred_block[0:16383],
    output reg                     valid,
    input  wire                    ready
);

// Simple stub: return constant pixels
reg busy;
reg [3:0] cycle_count;

integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        valid <= 1'b0;
        busy <= 1'b0;
        cycle_count <= 4'd0;
        ref_read_en <= 1'b0;
        ref_read_addr <= 32'd0;
        for (i = 0; i < 16384; i = i + 1) begin
            pred_block[i] <= 10'd128;
        end
    end else begin
        if (start && !busy) begin
            busy <= 1'b1;
            cycle_count <= 4'd0;
            ref_read_addr <= block_x;
            ref_read_en <= 1'b1;
        end else if (busy) begin
            if (cycle_count < 4'd4) begin
                cycle_count <= cycle_count + 1;
                ref_read_en <= 1'b0;
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