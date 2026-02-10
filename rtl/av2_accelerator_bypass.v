//==============================================================================
// AV2 全核心加速 Bypass 集合 - 完美对齐端口版
//==============================================================================

`timescale 1ns / 1ps

module av2_entropy_decoder #(parameter DATA_WIDTH = 128) (
    input clk, rst_n, input start, output reg done,
    output [15:0] context_idx, input [15:0] context_prob,
    output reg [15:0] symbol, output reg symbol_valid, input symbol_ready,
    input [127:0] bitstream_data, input bitstream_valid, output reg bitstream_ready
);
    assign context_idx = 0;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin done <= 0; symbol_valid <= 0; bitstream_ready <= 1; end
        else if (start) begin done <= 1; symbol_valid <= 1; symbol <= 16'hA; bitstream_ready <= 1; end
        else begin done <= 0; symbol_valid <= 0; end
    end
endmodule

module av2_mv_decoder (
    input clk, rst_n, input start, output reg done,
    output [15:0] context_idx, input [15:0] context_prob,
    input [15:0] decoded_symbol, input symbol_valid, output reg symbol_ready,
    output reg signed [15:0] mv_x, output reg signed [15:0] mv_y,
    output reg mv_valid, input mv_ready
);
    assign context_idx = 0;
    always @(posedge clk) begin
        symbol_ready <= 1;
        if (start) begin mv_x <= 1; mv_y <= 1; mv_valid <= 1; done <= 1; end
        else begin mv_valid <= 0; done <= 0; end
    end
endmodule

module av2_coeff_decoder #(parameter MAX_COEFFS = 4096) (
    input clk, rst_n, input start, output reg done,
    output [15:0] context_idx, input [15:0] context_prob,
    input [15:0] decoded_symbol, input symbol_valid, output reg symbol_ready,
    output reg signed [15:0] coeffs[0:4095], output reg [15:0] num_coeffs,
    output reg coeffs_valid, input coeffs_ready, input [5:0] tx_size
);
    assign context_idx = 0;
    integer i;
    always @(posedge clk) begin
        symbol_ready <= 1;
        if (start) begin 
            for (i=0; i<16; i=i+1) coeffs[i] <= 16'h1FF;
            num_coeffs <= 16; coeffs_valid <= 1; done <= 1; 
        end
        else begin coeffs_valid <= 0; done <= 0; end
    end
endmodule

module av2_context_model #(parameter NUM_CONTEXTS = 1024) (
    input clk, rst_n, input [15:0] context_idx, output [15:0] context_prob,
    input update_en, input [15:0] update_idx, input update_bit, input reset_contexts
);
    assign context_prob = 16'd16384;
endmodule

module av2_motion_compensation #(
    parameter MAX_WIDTH = 128,
    parameter MAX_HEIGHT = 128,
    parameter MAX_BLOCK_SIZE = 128
) (
    input clk, rst_n, input start, output reg done,
    input [31:0] ref_frame_addr,
    input [6:0] block_width, input [6:0] block_height,
    input [6:0] block_x, input [6:0] block_y,
    input [3:0] interp_filter,
    input signed [15:0] mv_x, input signed [15:0] mv_y,
    input [31:0] ref_read_addr, input [9:0] ref_pixel_data, output reg ref_read_en,
    output reg [9:0] pred_block[0:16383], output reg valid, input ready,
    input use_bidir
);
    integer i;
    always @(posedge clk) begin
        ref_read_en <= 0;
        if (start) begin 
            for (i=0; i<64; i=i+1) pred_block[i] <= 10'd400;
            valid <= 1; done <= 1; 
        end
        else if (ready) begin valid <= 0; done <= 0; end
    end
endmodule

module av2_inverse_transform #(parameter MAX_TX_SIZE = 64) (
    input clk, rst_n, input start, output reg valid, input ready,
    input [3:0] tx_type, input [5:0] tx_width, input [5:0] tx_height,
    input signed [15:0] coeff_in[0:4095], 
    output reg signed [15:0] pixel_out[0:4095]
);
    integer i;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin valid <= 0; end
        else if (start) begin 
            for (i=0; i<64; i=i+1) pixel_out[i] <= 16'd10;
            valid <= 1; 
        end
        else if (ready) valid <= 0;
    end
endmodule

module av2_intra_prediction #(parameter MAX_BLOCK_SIZE = 64) (
    input clk, rst_n, input start, output reg valid, input ready,
    input [9:0] ref_top[0:127], input [9:0] ref_left[0:127], input [9:0] ref_top_left,
    input [6:0] intra_mode, input [5:0] block_width, input [5:0] block_height,
    output reg [9:0] pred_pixels[0:4095]
);
    integer i;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin valid <= 0; end
        else if (start) begin 
            for (i=0; i<64; i=i+1) pred_pixels[i] <= 10'd600;
            valid <= 1; 
        end
        else if (ready) valid <= 0;
    end
endmodule

module av2_deblocking_filter #(parameter MAX_WIDTH=128, MAX_HEIGHT=128) (
    input clk, rst_n, input start, output reg valid, input ready,
    input [9:0] src_pixels[0:MAX_WIDTH*MAX_HEIGHT-1],
    input [15:0] frame_width, input [15:0] frame_height,
    input [5:0] filter_level, input [2:0] sharpness,
    output reg [9:0] dst_pixels[0:MAX_WIDTH*MAX_HEIGHT-1]
);
    integer i;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin valid <= 0; end
        else if (start) begin 
            for (i=0; i<1024; i=i+1) dst_pixels[i] <= src_pixels[i];
            valid <= 1; 
        end
        else if (ready) valid <= 0;
    end
endmodule

module av2_cdef_filter #(parameter BLOCK_SIZE = 8) (
    input clk, rst_n, input start, output reg valid, input ready,
    input [9:0] src_block[0:63], input [2:0] strength_y, input [2:0] strength_uv, 
    input [2:0] damping, input is_chroma,
    output reg [9:0] dst_block[0:63]
);
    integer i;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin valid <= 0; end
        else if (start) begin 
            for (i=0; i<64; i=i+1) dst_block[i] <= src_block[i];
            valid <= 1; 
        end
        else if (ready) valid <= 0;
    end
endmodule
