//==============================================================================
// Simple AV2 Decoder Testbench
// Quick verification of basic functionality
//==============================================================================

`timescale 1ns / 1ps

module tb_simple;

//==============================================================================
// Parameters
//==============================================================================

parameter CLK_PERIOD = 10;
parameter MAX_WIDTH  = 64;
parameter MAX_HEIGHT = 64;

//==============================================================================
// Signals
//==============================================================================

reg clk;
reg rst_n;

reg  [127:0] s_axis_tdata;
reg          s_axis_tvalid;
wire         s_axis_tready;
reg          s_axis_tlast;
reg  [15:0]  s_axis_tkeep;

wire [127:0] m_axis_tdata;
wire         m_axis_tvalid;
reg          m_axis_tready;
wire         m_axis_tlast;
wire [15:0]  m_axis_tkeep;
wire [1:0]   m_axis_tuser;

reg  [15:0]  reg_addr;
reg  [31:0]  reg_wdata;
reg          reg_wr_en;
reg          reg_rd_en;
wire [31:0]  reg_rdata;
wire         reg_ready;

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

wire         irq_frame_done;
wire         irq_error;

//==============================================================================
// DUT
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
// Clock
//==============================================================================

initial begin
    clk = 0;
    forever #(CLK_PERIOD/2) clk = ~clk;
end

//==============================================================================
// Simplified AXI Memory Model
//==============================================================================

reg [127:0] mem[0:1023];

always @(posedge clk) begin
    if (!rst_n) begin
        m_axi_awready <= 0;
        m_axi_wready <= 0;
        m_axi_bvalid <= 0;
        m_axi_arready <= 0;
        m_axi_rvalid <= 0;
    end else begin
        m_axi_awready <= m_axi_awvalid;
        m_axi_wready <= m_axi_wvalid;
        m_axi_arready <= m_axi_arvalid;
        
        if (m_axi_wvalid && m_axi_wready)
            mem[m_axi_awaddr[11:2]] <= m_axi_wdata;
        
        m_axi_bvalid <= m_axi_wvalid && m_axi_wready;
        m_axi_bresp <= 2'b00;
        
        if (m_axi_arvalid && m_axi_arready) begin
            m_axi_rdata <= mem[m_axi_araddr[11:2]];
            m_axi_rvalid <= 1;
            m_axi_rlast <= 1;
            m_axi_rresp <= 2'b00;
        end else if (m_axi_rvalid && m_axi_rready) begin
            m_axi_rvalid <= 0;
        end
    end
end

//==============================================================================
// Test
//==============================================================================

reg [31:0] rdata;
integer cycle_count;

initial begin
    // Initialize
    rst_n = 0;
    s_axis_tdata = 0;
    s_axis_tvalid = 0;
    s_axis_tlast = 0;
    s_axis_tkeep = 16'hFFFF;
    m_axis_tready = 1;
    reg_addr = 0;
    reg_wdata = 0;
    reg_wr_en = 0;
    reg_rd_en = 0;
    cycle_count = 0;
    
    // Reset
    #(CLK_PERIOD * 10);
    rst_n = 1;
    #(CLK_PERIOD * 5);
    
    $display("========================================");
    $display("Simple AV2 Decoder Test");
    $display("========================================");
    
    // Test 1: Read ID
    $display("\n[Test 1] Reading Decoder ID...");
    read_reg(16'h0000, rdata);
    $display("  ID: 0x%08X %s", rdata, 
        (rdata == 32'h41563200) ? "✓ PASS" : "✗ FAIL");
    
    // Test 2: Check status
    $display("\n[Test 2] Checking status...");
    read_reg(16'h0004, rdata);
    $display("  Busy: %d", rdata[0]);
    
    // Test 3: Send minimal data
    $display("\n[Test 3] Sending minimal data...");
    s_axis_tdata = {8'h0A, 32'h00000010, 88'h0};  // OBU header
    s_axis_tvalid = 1;
    s_axis_tlast = 1;
    @(posedge clk);
    while (!s_axis_tready) @(posedge clk);
    s_axis_tvalid = 0;
    s_axis_tlast = 0;
    $display("  Data sent");
    
    // Wait a bit
    #(CLK_PERIOD * 100);
    
    // Test 4: Check state
    $display("\n[Test 4] Checking decoder state...");
    read_reg(16'h0018, rdata);
    $display("  State: %d", rdata[3:0]);
    
    // Test 5: Read stats
    $display("\n[Test 5] Reading statistics...");
    read_reg(16'h0008, rdata);
    $display("  Frames decoded: %d", rdata);
    
    read_reg(16'h000C, rdata);
    $display("  Error count: %d", rdata);
    
    $display("\n========================================");
    $display("Test completed after %d cycles", cycle_count);
    $display("========================================");
    
    #(CLK_PERIOD * 10);
    $finish;
end

always @(posedge clk) begin
    if (rst_n)
        cycle_count <= cycle_count + 1;
end

task read_reg;
    input [15:0] addr;
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

endmodule