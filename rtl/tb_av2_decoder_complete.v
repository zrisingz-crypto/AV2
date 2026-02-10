//==============================================================================
// AV2 解码器完整测试平台
// 测试完整的解码流程
//==============================================================================

`timescale 1ns / 1ps

module tb_av2_decoder_complete;

//==============================================================================
// 参数定义
//==============================================================================

parameter CLK_PERIOD = 10;  // 100MHz
parameter MAX_WIDTH  = 128;
parameter MAX_HEIGHT = 128;

//==============================================================================
// 信号定义
//==============================================================================

// 时钟和复位
reg clk;
reg rst_n;

// AXI4-Stream 输入
reg  [127:0] s_axis_tdata;
reg          s_axis_tvalid;
wire         s_axis_tready;
reg          s_axis_tlast;
reg  [15:0]  s_axis_tkeep;

// AXI4-Stream 输出
wire [127:0] m_axis_tdata;
wire         m_axis_tvalid;
reg          m_axis_tready;
wire         m_axis_tlast;
wire [15:0]  m_axis_tkeep;
wire [1:0]   m_axis_tuser;

// 控制寄存器
reg  [15:0]  reg_addr;
reg  [31:0]  reg_wdata;
reg          reg_wr_en;
reg          reg_rd_en;
wire [31:0]  reg_rdata;
wire         reg_ready;

// AXI4 帧缓冲接口（简化为 BRAM）
wire [31:0]  m_axi_awaddr;
wire [7:0]   m_axi_awlen;
wire         m_axi_awvalid;
reg          m_axi_awready;
wire [127:0] m_axi_wdata;
wire         m_axi_wlast;
wire         m_axi_wvalid;
reg          m_axi_wready;
reg  [1:0]   m_axi_bresp;
reg          m_axi_bvalid;
wire         m_axi_bready;

wire [31:0]  m_axi_araddr;
wire [7:0]   m_axi_arlen;
wire         m_axi_arvalid;
reg          m_axi_arready;
reg  [127:0] m_axi_rdata;
reg          m_axi_rlast;
reg  [1:0]   m_axi_rresp;
reg          m_axi_rvalid;
wire         m_axi_rready;

// 中断
wire         irq_frame_done;
wire         irq_error;

//==============================================================================
// DUT 实例化
//==============================================================================

av2_decoder_top #(
    .MAX_FRAME_WIDTH (MAX_WIDTH),
    .MAX_FRAME_HEIGHT(MAX_HEIGHT),
    .PIXEL_WIDTH     (10),
    .AXI_DATA_WIDTH  (128),
    .AXI_ADDR_WIDTH  (32)
) dut (
    .clk             (clk),
    .rst_n           (rst_n),
    
    .s_axis_tdata    (s_axis_tdata),
    .s_axis_tvalid   (s_axis_tvalid),
    .s_axis_tready   (s_axis_tready),
    .s_axis_tlast    (s_axis_tlast),
    .s_axis_tkeep    (s_axis_tkeep),
    
    .m_axis_tdata    (m_axis_tdata),
    .m_axis_tvalid   (m_axis_tvalid),
    .m_axis_tready   (m_axis_tready),
    .m_axis_tlast    (m_axis_tlast),
    .m_axis_tkeep    (m_axis_tkeep),
    .m_axis_tuser    (m_axis_tuser),
    
    .reg_addr        (reg_addr),
    .reg_wdata       (reg_wdata),
    .reg_wr_en       (reg_wr_en),
    .reg_rd_en       (reg_rd_en),
    .reg_rdata       (reg_rdata),
    .reg_ready       (reg_ready),
    
    .m_axi_awaddr    (m_axi_awaddr),
    .m_axi_awlen     (m_axi_awlen),
    .m_axi_awvalid   (m_axi_awvalid),
    .m_axi_awready   (m_axi_awready),
    .m_axi_wdata     (m_axi_wdata),
    .m_axi_wlast     (m_axi_wlast),
    .m_axi_wvalid    (m_axi_wvalid),
    .m_axi_wready    (m_axi_wready),
    .m_axi_bresp     (m_axi_bresp),
    .m_axi_bvalid    (m_axi_bvalid),
    .m_axi_bready    (m_axi_bready),
    
    .m_axi_araddr    (m_axi_araddr),
    .m_axi_arlen     (m_axi_arlen),
    .m_axi_arvalid   (m_axi_arvalid),
    .m_axi_arready   (m_axi_arready),
    .m_axi_rdata     (m_axi_rdata),
    .m_axi_rlast     (m_axi_rlast),
    .m_axi_rresp     (m_axi_rresp),
    .m_axi_rvalid    (m_axi_rvalid),
    .m_axi_rready    (m_axi_rready),
    
    .irq_frame_done  (irq_frame_done),
    .irq_error       (irq_error)
);

//==============================================================================
// 时钟生成
//==============================================================================

initial begin
    clk = 0;
    forever #(CLK_PERIOD/2) clk = ~clk;
end

//==============================================================================
// 简化的帧缓冲模型（BRAM）
//==============================================================================

reg [127:0] frame_buffer[0:4095];  // 简化的帧缓冲

// AXI4 写响应
always @(posedge clk) begin
    if (!rst_n) begin
        m_axi_awready <= 1'b0;
        m_axi_wready <= 1'b0;
        m_axi_bvalid <= 1'b0;
    end else begin
        m_axi_awready <= m_axi_awvalid;
        m_axi_wready <= m_axi_wvalid;
        
        if (m_axi_wvalid && m_axi_wready) begin
            frame_buffer[m_axi_awaddr[13:2]] <= m_axi_wdata;
            m_axi_bvalid <= 1'b1;
            m_axi_bresp <= 2'b00;  // OKAY
        end else if (m_axi_bvalid && m_axi_bready) begin
            m_axi_bvalid <= 1'b0;
        end
    end
end

// AXI4 读响应
always @(posedge clk) begin
    if (!rst_n) begin
        m_axi_arready <= 1'b0;
        m_axi_rvalid <= 1'b0;
    end else begin
        m_axi_arready <= m_axi_arvalid;
        
        if (m_axi_arvalid && m_axi_arready) begin
            m_axi_rdata <= frame_buffer[m_axi_araddr[13:2]];
            m_axi_rvalid <= 1'b1;
            m_axi_rlast <= 1'b1;
            m_axi_rresp <= 2'b00;  // OKAY
        end else if (m_axi_rvalid && m_axi_rready) begin
            m_axi_rvalid <= 1'b0;
        end
    end
end

//==============================================================================
// 测试激励
//==============================================================================

integer i;
reg [31:0] test_data;

initial begin
    // 初始化
    rst_n = 0;
    s_axis_tdata = 128'h0;
    s_axis_tvalid = 0;
    s_axis_tlast = 0;
    s_axis_tkeep = 16'hFFFF;
    m_axis_tready = 1;
    reg_addr = 16'h0;
    reg_wdata = 32'h0;
    reg_wr_en = 0;
    reg_rd_en = 0;
    
    // 复位
    #(CLK_PERIOD * 10);
    rst_n = 1;
    #(CLK_PERIOD * 5);
    
    $display("========================================");
    $display("AV2 Decoder Test Started");
    $display("========================================");
    
    // 测试 1: 读取解码器 ID
    $display("\n[Test 1] Reading Decoder ID...");
    read_register(16'h0000, test_data);
    if (test_data == 32'h41563200) begin
        $display("✓ Decoder ID correct: 0x%08X (AV2)", test_data);
    end else begin
        $display("✗ Decoder ID incorrect: 0x%08X", test_data);
    end
    
    // 测试 2: 检查解码器状态
    $display("\n[Test 2] Checking Decoder Status...");
    read_register(16'h0004, test_data);
    $display("  Decoder Busy: %d", test_data[0]);
    
    // 测试 3: 发送简单的 OBU 数据
    $display("\n[Test 3] Sending Test Bitstream...");
    
    // 发送 OBU Sequence Header
    send_obu_packet(8'h0A, 32'h00000010, 0);  // OBU type=1 (SEQUENCE_HEADER), size=16
    
    // 发送 OBU Frame Header  
    send_obu_packet(8'h1A, 32'h00000020, 0);  // OBU type=3 (FRAME_HEADER), size=32
    
    // 发送帧数据
    for (i = 0; i < 10; i = i + 1) begin
        s_axis_tdata = {i[7:0], 120'h0};
        s_axis_tvalid = 1;
        s_axis_tlast = (i == 9);
        @(posedge clk);
        while (!s_axis_tready) @(posedge clk);
    end
    s_axis_tvalid = 0;
    s_axis_tlast = 0;
    
    $display("  Bitstream sent");
    
    // 等待解码完成
    $display("\n[Test 4] Waiting for Frame Decode...");
    wait(irq_frame_done);
    $display("  ✓ Frame decode interrupt received!");
    
    // 读取统计信息
    $display("\n[Test 5] Reading Statistics...");
    read_register(16'h0008, test_data);
    $display("  Frames Decoded: %d", test_data);
    
    read_register(16'h000C, test_data);
    $display("  Error Count: %d", test_data);
    
    read_register(16'h0010, test_data);
    $display("  Frame Width: %d", test_data[15:0]);
    
    read_register(16'h0014, test_data);
    $display("  Frame Height: %d", test_data[15:0]);
    
    // 读取输出帧
    $display("\n[Test 6] Reading Decoded Frame...");
    i = 0;
    while (!m_axis_tlast && i < 1000) begin
        if (m_axis_tvalid && m_axis_tready) begin
            if (i < 5) begin
                $display("  Output[%0d]: 0x%032X", i, m_axis_tdata);
            end
            i = i + 1;
        end
        @(posedge clk);
    end
    $display("  Total output packets: %d", i);
    
    // 测试完成
    #(CLK_PERIOD * 100);
    
    $display("\n========================================");
    $display("AV2 Decoder Test Completed");
    $display("========================================");
    
    $finish;
end

//==============================================================================
// 任务定义
//==============================================================================

// 读寄存器
task read_register;
    input  [15:0] addr;
    output [31:0] data;
    begin
        @(posedge clk);
        reg_addr = addr;
        reg_rd_en = 1;
        @(posedge clk);
        while (!reg_ready) @(posedge clk);
        data = reg_rdata;
        reg_rd_en = 0;
        @(posedge clk);
    end
endtask

// 写寄存器
task write_register;
    input [15:0] addr;
    input [31:0] data;
    begin
        @(posedge clk);
        reg_addr = addr;
        reg_wdata = data;
        reg_wr_en = 1;
        @(posedge clk);
        while (!reg_ready) @(posedge clk);
        reg_wr_en = 0;
        @(posedge clk);
    end
endtask

// 发送 OBU 包
task send_obu_packet;
    input [7:0]  obu_header;
    input [31:0] obu_size;
    input        is_last;
    begin
        @(posedge clk);
        s_axis_tdata = {obu_header, obu_size[7:0], 112'h0};
        s_axis_tvalid = 1;
        s_axis_tlast = is_last;
        @(posedge clk);
        while (!s_axis_tready) @(posedge clk);
        s_axis_tvalid = 0;
        s_axis_tlast = 0;
    end
endtask

//==============================================================================
// 波形转储
//==============================================================================

initial begin
    $dumpfile("av2_decoder.vcd");
    $dumpvars(0, tb_av2_decoder_complete);
end

//==============================================================================
// 超时保护
//==============================================================================

initial begin
    #(CLK_PERIOD * 100000);  // 1ms 超时
    $display("\n========================================");
    $display("ERROR: Simulation Timeout!");
    $display("========================================");
    $finish;
end

endmodule
