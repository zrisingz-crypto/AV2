//==============================================================================
// AV2 环路滤波器集合 - 修复 Icarus Verilog 兼容性版本
//==============================================================================

`timescale 1ns / 1ps

//==============================================================================
// Deblocking Filter
//==============================================================================

module av2_deblocking_filter #(
    parameter MAX_WIDTH = 128,
    parameter MAX_HEIGHT = 128
)(
    input  wire        clk,
    input  wire        rst_n,
    input  wire [9:0]  src_pixels[0:MAX_WIDTH*MAX_HEIGHT-1],
    input  wire [15:0] frame_width,
    input  wire [15:0] frame_height,
    input  wire [5:0]  filter_level,
    input  wire [2:0]  sharpness,
    input  wire        start,
    output reg [9:0]   dst_pixels[0:MAX_WIDTH*MAX_HEIGHT-1],
    output reg         valid,
    input  wire        ready
);

localparam IDLE = 2'd0;
localparam FILTER = 2'd1;
localparam OUTPUT = 2'd2;

reg [1:0] state;
reg [15:0] row, col;
integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        valid <= 1'b0;
    end else begin
        case (state)
            IDLE: begin
                valid <= 1'b0;
                if (start) begin
                    for (i = 0; i < MAX_WIDTH*MAX_HEIGHT; i = i + 1) dst_pixels[i] <= src_pixels[i];
                    state <= FILTER;
                end
            end
            FILTER: begin
                // 目前仅做直通
                valid <= 1'b1;
                state <= OUTPUT;
            end
            OUTPUT: begin
                if (ready) begin
                    valid <= 1'b0;
                    state <= IDLE;
                end
            end
        endcase
    end
end
endmodule

//==============================================================================
// CDEF Filter
//==============================================================================

module av2_cdef_filter #(
    parameter BLOCK_SIZE = 8
)(
    input  wire        clk,
    input  wire        rst_n,
    input  wire [9:0]  src_block[0:BLOCK_SIZE*BLOCK_SIZE-1],
    input  wire [2:0]  strength_y,
    input  wire [2:0]  strength_uv,
    input  wire [2:0]  damping,
    input  wire        is_chroma,
    input  wire        start,
    output reg [9:0]   dst_block[0:BLOCK_SIZE*BLOCK_SIZE-1],
    output reg         valid,
    input  wire        ready
);

reg [1:0] state;
reg [3:0] r, c;
reg signed [15:0] diff;
integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= 2'd0;
        valid <= 1'b0;
    end else begin
        case (state)
            2'd0: begin
                valid <= 1'b0;
                if (start) begin
                    r <= 4'd0;
                    c <= 4'd0;
                    state <= 2'd1;
                end
            end
            2'd1: begin
                if (r < BLOCK_SIZE) begin
                    if (c < BLOCK_SIZE) begin
                        // 简化 CDEF: 仅做简单的平滑
                        diff = 0;
                        if (c > 0) diff = diff + src_block[r * BLOCK_SIZE + c - 1] - src_block[r * BLOCK_SIZE + c];
                        if (c < 7) diff = diff + src_block[r * BLOCK_SIZE + c + 1] - src_block[r * BLOCK_SIZE + c];
                        dst_block[r * BLOCK_SIZE + c] <= src_block[r * BLOCK_SIZE + c] + (diff >> 3);
                        c <= c + 1;
                    end else begin
                        c <= 4'd0;
                        r <= r + 1;
                    end
                end else begin
                    valid <= 1'b1;
                    state <= 2'd2;
                end
            end
            2'd2: begin
                if (ready) begin
                    valid <= 1'b0;
                    state <= 2'd0;
                end
            end
        endcase
    end
end
endmodule

//==============================================================================
// Loop Restoration
//==============================================================================

module av2_loop_restoration #(
    parameter BLOCK_SIZE = 64
)(
    input  wire        clk,
    input  wire        rst_n,
    input  wire [9:0]  src_block[0:BLOCK_SIZE*BLOCK_SIZE-1],
    input  wire signed [7:0] wiener_coeffs[0:48],
    input  wire        start,
    output reg [9:0]   dst_block[0:BLOCK_SIZE*BLOCK_SIZE-1],
    output reg         valid,
    input  wire        ready
);

reg [1:0] state;
integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= 2'd0;
        valid <= 1'b0;
    end else begin
        case (state)
            2'd0: begin
                valid <= 1'b0;
                if (start) begin
                    for (i = 0; i < BLOCK_SIZE*BLOCK_SIZE; i = i + 1) dst_block[i] <= src_block[i];
                    state <= 2'd1;
                end
            end
            2'd1: begin
                valid <= 1'b1;
                state <= 2'd2;
            end
            2'd2: begin
                if (ready) begin
                    valid <= 1'b0;
                    state <= 2'd0;
                end
            end
        endcase
    end
end
endmodule
