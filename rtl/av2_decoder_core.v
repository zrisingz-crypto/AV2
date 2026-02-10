//==============================================================================
// AV2 OBU 解析器 - 修复位偏移 (数据在 AXI 总线高位)
//==============================================================================

`timescale 1ns / 1ps

module av2_obu_parser #(
    parameter DATA_WIDTH = 128
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire [DATA_WIDTH-1:0]   s_axis_tdata,
    input  wire                     s_axis_tvalid,
    output wire                     s_axis_tready,
    input  wire                     s_axis_tlast,
    output reg  [3:0]               obu_type,
    output reg  [31:0]              obu_size,
    output reg                      obu_valid,
    input  wire                     obu_ready
);

localparam IDLE = 2'd0;
localparam DONE = 2'd1;

reg [1:0] state;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        obu_type <= 4'h0;
        obu_size <= 32'h0;
        obu_valid <= 1'b0;
    end else begin
        case (state)
            IDLE: begin
                if (s_axis_tvalid) begin
                    // 改回从低位解析 (对应单字节推送)
                    obu_type <= s_axis_tdata[6:3];
                    obu_size <= {24'h0, s_axis_tdata[15:8]}; 
                    obu_valid <= 1'b1;
                    state <= DONE;
                end
            end
            
            DONE: begin
                if (obu_ready) begin
                    obu_valid <= 1'b0;
                    state <= IDLE;
                end
            end
        endcase
    end
end

assign s_axis_tready = (state == IDLE);

endmodule


//==============================================================================
// AV2 帧头解析器 - 修复位偏移
//==============================================================================
module av2_frame_header_parser (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        obu_valid,
    input  wire [127:0] obu_data,
    output reg  [1:0]  frame_type,
    output reg  [15:0] frame_width,
    output reg  [15:0] frame_height,
    output reg  [7:0]  qindex,
    output reg         header_valid,
    input  wire        header_ready
);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            header_valid <= 0;
            frame_type <= 2'b0;  // 默认为 KEY_FRAME
            frame_width <= 64;
            frame_height <= 64;
        end else if (obu_valid) begin
$display("    [DEBUG] Frame Header Raw: %h", obu_data[63:0]);
            // 数据紧跟在 OBU Header (1B) 和 Size (1B) 之后，即从 16 位开始
            frame_type   <= obu_data[17:16];
            frame_width  <= {8'h0, obu_data[31:24]};
            frame_height <= {8'h0, obu_data[39:32]};
            qindex       <= 8'd128;
            header_valid <= 1'b1;
        end else if (header_ready) begin
            header_valid <= 1'b0;
        end
    end
endmodule

// 重建与滤波占位模块保持不变...
module av2_reconstruction (
    input clk, rst_n, start, input [1:0] frame_type,
    output reg [31:0] fb_rd_addr, input [127:0] fb_rd_data, output reg fb_rd_en,
    output reg recon_valid, input recon_ready, output reg recon_done
);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin recon_valid <= 0; recon_done <= 0; end
        else if (start) begin recon_valid <= 1; if (recon_ready) recon_done <= 1; end
        else begin recon_valid <= 0; recon_done <= 0; end
    end
endmodule

module av2_loop_filter (
    input clk, rst_n, start, input [15:0] frame_width, frame_height,
    output reg [31:0] fb_addr, input [127:0] fb_data,
    output reg filter_valid, input filter_ready, output reg filter_done
);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin filter_valid <= 0; filter_done <= 0; end
        else if (start) begin filter_valid <= 1; if (filter_ready) filter_done <= 1; end
        else begin filter_valid <= 0; filter_done <= 0; end
    end
endmodule
