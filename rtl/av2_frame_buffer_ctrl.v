//==============================================================================
// Frame Buffer Controller Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_frame_buffer_ctrl #(
    parameter ADDR_WIDTH = 32,
    parameter DATA_WIDTH = 128
)(
    input  wire                    clk,
    input  wire                    rst_n,
    
    // Internal Interface
    input  wire [ADDR_WIDTH-1:0]  wr_addr,
    input  wire [DATA_WIDTH-1:0]  wr_data,
    input  wire                    wr_en,
    input  wire [ADDR_WIDTH-1:0]  rd_addr,
    output reg  [DATA_WIDTH-1:0]  rd_data,
    input  wire                    rd_en,
    input  wire [1:0]              rd_sel_plane,  // 0: Y, 1: U, 2: V
    
    // AXI4 Master Interface (simplified)
    output wire [ADDR_WIDTH-1:0]  m_axi_awaddr,
    output wire [7:0]              m_axi_awlen,
    output wire                    m_axi_awvalid,
    input  wire                    m_axi_awready,
    output wire [DATA_WIDTH-1:0]   m_axi_wdata,
    output wire                    m_axi_wlast,
    output wire                    m_axi_wvalid,
    input  wire                    m_axi_wready,
    input  wire [1:0]              m_axi_bresp,
    input  wire                    m_axi_bvalid,
    output wire                    m_axi_bready,
    
    output wire [ADDR_WIDTH-1:0]  m_axi_araddr,
    output wire [7:0]              m_axi_arlen,
    output wire                    m_axi_arvalid,
    input  wire                    m_axi_arready,
    input  wire [DATA_WIDTH-1:0]   m_axi_rdata,
    input  wire                    m_axi_rlast,
    input  wire [1:0]              m_axi_rresp,
    input  wire                    m_axi_rvalid,
    output wire                    m_axi_rready
);

// Internal frame buffer - enough for 64x64 frame (4096 pixels) in YUV420 format
// Y: 4096 bytes = 256 words (word addresses 0-255)
// U: 1024 bytes = 64 words (word addresses 256-319)
// V: 1024 bytes = 64 words (word addresses 320-383)
// Total: 384 words
reg [DATA_WIDTH-1:0] frame_buffer[0:383];

// Write address is already in word format from tile decoder (0-383)
// Tile decoder provides absolute word addresses, no conversion needed

// Plane offset calculation for YUV420 (in word addresses)
wire [8:0] y_offset = 9'd0;
wire [8:0] u_offset = 9'd256;   // U plane starts at word 256
wire [8:0] v_offset = 9'd320;   // V plane starts at word 320

wire [8:0] read_addr_offset = 
    (rd_sel_plane == 2'd0) ? y_offset :
    (rd_sel_plane == 2'd1) ? u_offset :
    (rd_sel_plane == 2'd2) ? v_offset : y_offset;

reg [8:0] actual_read_addr;

// Write operation
always @(posedge clk) begin
    if (wr_en) begin
        // wr_addr is already in word format, use directly
        frame_buffer[wr_addr[8:0]] <= wr_data;
    end
end

// Read operation (corrected timing)
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        rd_data <= {DATA_WIDTH{1'b0}};
        actual_read_addr <= 9'd0;
    end else if (rd_en) begin
        // Read address is relative to current plane, need to apply offset
        actual_read_addr <= rd_addr[8:0] + read_addr_offset;
        rd_data <= frame_buffer[rd_addr[8:0] + read_addr_offset];  // Use new address immediately
    end
end

// AXI4 interface stubs (connect to testbench)
assign m_axi_awaddr = wr_addr;
assign m_axi_awlen = 8'd0;
assign m_axi_awvalid = wr_en;
assign m_axi_wdata = wr_data;
assign m_axi_wlast = 1'b1;
assign m_axi_wvalid = wr_en;
assign m_axi_bready = 1'b1;

assign m_axi_araddr = rd_addr;
assign m_axi_arlen = 8'd0;
assign m_axi_arvalid = rd_en;
assign m_axi_rready = 1'b1;

endmodule