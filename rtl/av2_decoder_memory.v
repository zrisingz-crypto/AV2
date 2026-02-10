//==============================================================================
// AV2 Frame Buffer Controller with AXI4 Interface
//==============================================================================

`timescale 1ns / 1ps

module av2_frame_buffer_ctrl #(
    parameter ADDR_WIDTH = 32,
    parameter DATA_WIDTH = 128
)(
    input  wire                     clk,
    input  wire                     rst_n,
    
    // Internal Write Interface
    input  wire [ADDR_WIDTH-1:0]   wr_addr,
    input  wire [DATA_WIDTH-1:0]   wr_data,
    input  wire                     wr_en,
    
    // Internal Read Interface
    input  wire [ADDR_WIDTH-1:0]   rd_addr,
    output reg  [DATA_WIDTH-1:0]   rd_data,
    input  wire                     rd_en,
    
    // AXI4 Master Write Interface
    output reg  [ADDR_WIDTH-1:0]   m_axi_awaddr,
    output reg  [7:0]               m_axi_awlen,
    output reg                      m_axi_awvalid,
    input  wire                     m_axi_awready,
    output reg  [DATA_WIDTH-1:0]   m_axi_wdata,
    output reg                      m_axi_wlast,
    output reg                      m_axi_wvalid,
    input  wire                     m_axi_wready,
    input  wire [1:0]               m_axi_bresp,
    input  wire                     m_axi_bvalid,
    output reg                      m_axi_bready,
    
    // AXI4 Master Read Interface
    output reg  [ADDR_WIDTH-1:0]   m_axi_araddr,
    output reg  [7:0]               m_axi_arlen,
    output reg                      m_axi_arvalid,
    input  wire                     m_axi_arready,
    input  wire [DATA_WIDTH-1:0]   m_axi_rdata,
    input  wire                     m_axi_rlast,
    input  wire [1:0]               m_axi_rresp,
    input  wire                     m_axi_rvalid,
    output reg                      m_axi_rready
);

// Write State Machine
localparam WR_IDLE  = 2'd0;
localparam WR_ADDR  = 2'd1;
localparam WR_DATA  = 2'd2;
localparam WR_RESP  = 2'd3;

reg [1:0] wr_state;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        wr_state <= WR_IDLE;
        m_axi_awaddr <= {ADDR_WIDTH{1'b0}};
        m_axi_awlen <= 8'd0;
        m_axi_awvalid <= 1'b0;
        m_axi_wdata <= {DATA_WIDTH{1'b0}};
        m_axi_wlast <= 1'b0;
        m_axi_wvalid <= 1'b0;
        m_axi_bready <= 1'b0;
    end else begin
        case (wr_state)
            WR_IDLE: begin
                if (wr_en) begin
                    m_axi_awaddr <= wr_addr;
                    m_axi_awlen <= 8'd0;  // Single transfer
                    m_axi_awvalid <= 1'b1;
                    wr_state <= WR_ADDR;
                end
            end
            
            WR_ADDR: begin
                if (m_axi_awready) begin
                    m_axi_awvalid <= 1'b0;
                    m_axi_wdata <= wr_data;
                    m_axi_wvalid <= 1'b1;
                    m_axi_wlast <= 1'b1;
                    wr_state <= WR_DATA;
                end
            end
            
            WR_DATA: begin
                if (m_axi_wready) begin
                    m_axi_wvalid <= 1'b0;
                    m_axi_wlast <= 1'b0;
                    m_axi_bready <= 1'b1;
                    wr_state <= WR_RESP;
                end
            end
            
            WR_RESP: begin
                if (m_axi_bvalid) begin
                    m_axi_bready <= 1'b0;
                    wr_state <= WR_IDLE;
                end
            end
        endcase
    end
end

// Read State Machine
localparam RD_IDLE = 2'd0;
localparam RD_ADDR = 2'd1;
localparam RD_DATA = 2'd2;

reg [1:0] rd_state;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        rd_state <= RD_IDLE;
        m_axi_araddr <= {ADDR_WIDTH{1'b0}};
        m_axi_arlen <= 8'd0;
        m_axi_arvalid <= 1'b0;
        m_axi_rready <= 1'b0;
        rd_data <= {DATA_WIDTH{1'b0}};
    end else begin
        case (rd_state)
            RD_IDLE: begin
                if (rd_en) begin
                    m_axi_araddr <= rd_addr;
                    m_axi_arlen <= 8'd0;  // Single transfer
                    m_axi_arvalid <= 1'b1;
                    rd_state <= RD_ADDR;
                end
            end
            
            RD_ADDR: begin
                if (m_axi_arready) begin
                    m_axi_arvalid <= 1'b0;
                    m_axi_rready <= 1'b1;
                    rd_state <= RD_DATA;
                end
            end
            
            RD_DATA: begin
                if (m_axi_rvalid) begin
                    rd_data <= m_axi_rdata;
                    m_axi_rready <= 1'b0;
                    rd_state <= RD_IDLE;
                end
            end
        endcase
    end
end

endmodule


//==============================================================================
// AV2 Output Controller
//==============================================================================

module av2_output_ctrl #(
    parameter DATA_WIDTH = 128
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     start,
    
    input  wire [15:0]              frame_width,
    input  wire [15:0]              frame_height,
    
    output reg  [31:0]              fb_rd_addr,
    input  wire [DATA_WIDTH-1:0]   fb_rd_data,
    
    output reg  [DATA_WIDTH-1:0]   m_axis_tdata,
    output reg                      m_axis_tvalid,
    input  wire                     m_axis_tready,
    output reg                      m_axis_tlast,
    output reg  [DATA_WIDTH/8-1:0]  m_axis_tkeep,
    output reg  [1:0]               m_axis_tuser
);

// Output state
localparam IDLE     = 2'd0;
localparam READING  = 2'd1;
localparam SENDING  = 2'd2;
localparam DONE     = 2'd3;

reg [1:0] state;
reg [31:0] pixel_count;
reg [31:0] total_pixels;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        fb_rd_addr <= 32'd0;
        m_axis_tdata <= {DATA_WIDTH{1'b0}};
        m_axis_tvalid <= 1'b0;
        m_axis_tlast <= 1'b0;
        m_axis_tkeep <= {(DATA_WIDTH/8){1'b1}};
        m_axis_tuser <= 2'b00;
        pixel_count <= 32'd0;
        total_pixels <= 32'd0;
    end else begin
        case (state)
            IDLE: begin
                if (start) begin
                    pixel_count <= 32'd0;
                    total_pixels <= frame_width * frame_height;
                    fb_rd_addr <= 32'd0;
                    m_axis_tuser <= 2'b10;  // Frame start
                    state <= READING;
                end
            end
            
            READING: begin
                // Read from frame buffer
                fb_rd_addr <= fb_rd_addr + (DATA_WIDTH / 8);
                m_axis_tdata <= fb_rd_data;
                m_axis_tvalid <= 1'b1;
                state <= SENDING;
            end
            
            SENDING: begin
                if (m_axis_tready) begin
                    pixel_count <= pixel_count + (DATA_WIDTH / 8);
                    
                    if (pixel_count >= total_pixels - (DATA_WIDTH / 8)) begin
                        m_axis_tlast <= 1'b1;
                        m_axis_tuser <= 2'b01;  // Frame end
                        state <= DONE;
                    end else begin
                        m_axis_tuser <= 2'b00;
                        state <= READING;
                    end
                end
            end
            
            DONE: begin
                if (m_axis_tready) begin
                    m_axis_tvalid <= 1'b0;
                    m_axis_tlast <= 1'b0;
                    m_axis_tuser <= 2'b00;
                    state <= IDLE;
                end
            end
        endcase
    end
end

endmodule
