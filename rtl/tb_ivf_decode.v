//==============================================================================
// IVF Bitstream Decode Testbench
// Decodes real IVF file and outputs YUV data
//==============================================================================

`timescale 1ns / 1ps

module tb_ivf_decode;

//==============================================================================
// Parameters
//==============================================================================

parameter CLK_PERIOD = 10;
parameter MAX_WIDTH  = 64;
parameter MAX_HEIGHT = 64;
parameter FRAME_WIDTH = 64;
parameter FRAME_HEIGHT = 64;
parameter FRAME_SIZE = FRAME_WIDTH * FRAME_HEIGHT;  // 4096 pixels

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
// Include IVF bitstream data
//==============================================================================

`include "ivf_bitstream_data.v"

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
// Simple AXI Memory Model for frame buffer
//==============================================================================

reg [7:0] frame_buffer [0:FRAME_SIZE*3/2-1];  // YUV420: Y + UV/4 (8-bit bytes)

always @(posedge clk) begin
    if (!rst_n) begin
        m_axi_awready <= 0;
        m_axi_wready <= 0;
        m_axi_bvalid <= 0;
        m_axi_arready <= 0;
    end else begin
        m_axi_awready <= m_axi_awvalid;
        m_axi_wready <= m_axi_wvalid;
        m_axi_arready <= m_axi_arvalid;
        
        // Handle write - store 16 bytes per 128-bit transfer
        if (m_axi_wvalid && m_axi_wready) begin
            integer i;
            integer addr;
            
            for (i = 0; i < 16; i = i + 1) begin
                addr = m_axi_awaddr[11:0] + i;
                if (addr < FRAME_SIZE*3/2) begin
                    frame_buffer[addr] <= m_axi_wdata[i*8 +: 8];
                end
            end
        end
        
        m_axi_bvalid <= m_axi_wvalid && m_axi_wready;
        m_axi_bresp <= 2'b00;
        
        // Handle read - read 16 bytes per 128-bit transfer
        if (m_axi_arvalid && m_axi_arready) begin
            integer i;
            m_axi_rvalid <= 1;
            m_axi_rlast <= 1;
            m_axi_rresp <= 2'b00;
            
            for (i = 0; i < 16; i = i + 1) begin
                m_axi_rdata[i*8 +: 8] <= frame_buffer[m_axi_araddr[11:0] + i];
            end
        end else if (m_axi_rvalid && m_axi_rready) begin
            m_axi_rvalid <= 0;
        end
    end
end

//==============================================================================
// YUV Output Collection
//==============================================================================

integer yuv_file;
integer pixel_count;
integer frames_decoded;
integer y_pixel_count;  // Counter for Y plane pixels
integer u_pixel_count;  // Counter for U plane pixels
integer v_pixel_count;  // Counter for V plane pixels
reg output_complete;    // Flag indicating output is complete

initial begin
    yuv_file = $fopen("out/decoded_output.yuv", "wb");
    pixel_count = 0;
    frames_decoded = 0;
    y_pixel_count = 0;
    u_pixel_count = 0;
    v_pixel_count = 0;
    output_complete = 0;
end

// Collect decoded pixels from output
always @(posedge clk) begin
    integer i;
    
    if (!rst_n) begin
        y_pixel_count = 0;
        u_pixel_count = 0;
        v_pixel_count = 0;
        pixel_count = 0;
    end else if (m_axis_tvalid && m_axis_tready) begin
        // Write pixels on every valid transfer
        for (i = 0; i < 16; i = i + 1) begin
            if (m_axis_tkeep[i]) begin
                $fwrite(yuv_file, "%c", m_axis_tdata[i*8+:8]);
                pixel_count = pixel_count + 1;
            end
        end
        
        // Check for frame end
        if (m_axis_tlast && m_axis_tuser[0]) begin
            frames_decoded = frames_decoded + 1;
            output_complete = 1;  // Set flag when output is complete
            $display("Frame %d decoded: %d pixels written", frames_decoded - 1, pixel_count);
            $display("Output complete at time %0t", $time);
            pixel_count = 0;
        end
    end
end

//==============================================================================
// Bitstream Sender
//==============================================================================

integer bitstream_ptr;
integer current_frame;

task send_bitstream;
    input [31:0] frame_idx;
    input [31:0] frame_size;
    
    reg [31:0] byte_ptr;
    reg [127:0] word;
    reg [15:0] keep;
    integer i;
    
    begin
        $display("Sending Frame %d: %d bytes", frame_idx, frame_size);
        
        byte_ptr = 0;
        
        while (byte_ptr < frame_size) begin
            // Pack 16 bytes per 128-bit word
            word = 0;
            keep = 0;
            
            for (i = 0; i < 16 && byte_ptr < frame_size; i = i + 1) begin
                if (frame_idx == 0)
                    word[i*8 +: 8] = frame_0_data[byte_ptr];
                else
                    word[i*8 +: 8] = frame_1_data[byte_ptr];
                keep[i] = 1'b1;
                byte_ptr = byte_ptr + 1;
            end
            
            // Send word
            s_axis_tdata = word;
            s_axis_tkeep = keep;
            s_axis_tvalid = 1;
            
            @(posedge clk);
            while (!s_axis_tready) @(posedge clk);
        end
        
        // Send last
        s_axis_tlast = 1;
        @(posedge clk);
        while (!s_axis_tready) @(posedge clk);
        s_axis_tvalid = 0;
        s_axis_tlast = 0;
        
        $display("  Frame %d sent complete", frame_idx);
    end
endtask

//==============================================================================
// Test
//==============================================================================

reg [31:0] rdata;

initial begin
    integer i;
    
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
    bitstream_ptr = 0;
    current_frame = 0;
    
    // Reset
    #(CLK_PERIOD * 10);
    rst_n = 1;
    #(CLK_PERIOD * 5);
    
    $display("========================================");
    $display("IVF Bitstream Decode Test");
    $display("========================================");
    
    // Test 1: Read ID
    $display("\n[Test 1] Reading Decoder ID...");
    read_reg(16'h0000, rdata);
    $display("  ID: 0x%08X", rdata);
    if (rdata == 32'h41563200)
        $display("  PASS");
    else
        $display("  FAIL");
    
    // Decode first frame for comparison with SW decoder
    // Frame 0 has the actual frame data (3392 bytes)
    // Frame 1 is just a small header (47 bytes) with no frame data
    $display("\n[Frame 0] Starting decode...");
    
    // Send first frame bitstream (real frame data)
    send_bitstream(0, FRAME_0_SIZE);
    
    // Wait for output to complete (wait for tlast with frame_end)
    $display("Waiting for output to complete...");
    wait(output_complete == 1);
    $display("Output completed at time %0t", $time);
    
    // Wait a bit longer to ensure all data is written
    #(CLK_PERIOD * 100);
    
    // Check status
    read_reg(16'h0004, rdata);
    $display("  Busy: %d", rdata[0]);
    
    read_reg(16'h0008, rdata);
    $display("  Frames decoded: %d", rdata);
    
    // Final status
    $display("\n[Final] Checking decoder status...");
    read_reg(16'h0008, rdata);
    $display("  Total frames decoded: %d", rdata);
    
    read_reg(16'h000C, rdata);
    $display("  Error count: %d", rdata);
    
    $display("\n========================================");
    $display("Test Complete");
    $display("========================================");
    
    #(CLK_PERIOD * 100);
    $fclose(yuv_file);
    $finish;
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