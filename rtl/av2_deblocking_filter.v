//==============================================================================
// Deblocking Filter Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_deblocking_filter #(
    parameter MAX_WIDTH = 128,
    parameter MAX_HEIGHT = 128
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [9:0]             src_pixels[0:MAX_WIDTH*MAX_HEIGHT-1],
    input  wire [15:0]            frame_width,
    input  wire [15:0]            frame_height,
    input  wire [5:0]             filter_level,
    input  wire [2:0]             sharpness,
    input  wire                    start,
    output reg  [9:0]             dst_pixels[0:MAX_WIDTH*MAX_HEIGHT-1],
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
        for (i = 0; i < MAX_WIDTH * MAX_HEIGHT; i = i + 1) begin
            dst_pixels[i] <= 10'd0;
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