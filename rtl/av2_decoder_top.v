//==============================================================================
// AV2 Video Decoder - Top Level Module
// 
// Based on: Alliance for Open Media AVM Decoder
// Hardware Implementation: Simplified RTL Version
//
// Copyright (c) 2026
// License: BSD 3-Clause Clear License
//==============================================================================

`timescale 1ns / 1ps

module av2_decoder_top #(
    parameter MAX_FRAME_WIDTH  = 64,    // Reduced to 1 SB
    parameter MAX_FRAME_HEIGHT = 64,    // Reduced to 1 SB
    parameter PIXEL_WIDTH      = 10,    // 10-bit pixel depth
    parameter AXI_DATA_WIDTH   = 128,   // AXI bus width
    parameter AXI_ADDR_WIDTH   = 32
)(
    // Clock and Reset
    input  wire                         clk,
    input  wire                         rst_n,
    
    // AXI4-Stream Input Interface (Compressed Bitstream)
    input  wire [AXI_DATA_WIDTH-1:0]   s_axis_tdata,
    input  wire                         s_axis_tvalid,
    output wire                         s_axis_tready,
    input  wire                         s_axis_tlast,
    input  wire [AXI_DATA_WIDTH/8-1:0]  s_axis_tkeep,
    
    // AXI4-Stream Output Interface (Decoded Frame)
    output wire [AXI_DATA_WIDTH-1:0]   m_axis_tdata,
    output wire                         m_axis_tvalid,
    input  wire                         m_axis_tready,
    output wire                         m_axis_tlast,
    output wire [AXI_DATA_WIDTH/8-1:0]  m_axis_tkeep,
    output wire [1:0]                   m_axis_tuser,  // [1]: frame_start, [0]: frame_end
    
    // Control/Status Registers (APB/AXI-Lite Interface)
    input  wire [15:0]                  reg_addr,
    input  wire [31:0]                  reg_wdata,
    input  wire                         reg_wr_en,
    input  wire                         reg_rd_en,
    output wire [31:0]                  reg_rdata,
    output wire                         reg_ready,
    
    // Frame Buffer Memory Interface (AXI4 Master)
    output wire [AXI_ADDR_WIDTH-1:0]   m_axi_awaddr,
    output wire [7:0]                   m_axi_awlen,
    output wire                         m_axi_awvalid,
    input  wire                         m_axi_awready,
    output wire [AXI_DATA_WIDTH-1:0]   m_axi_wdata,
    output wire                         m_axi_wlast,
    output wire                         m_axi_wvalid,
    input  wire                         m_axi_wready,
    input  wire [1:0]                   m_axi_bresp,
    input  wire                         m_axi_bvalid,
    output wire                         m_axi_bready,
    
    output wire [AXI_ADDR_WIDTH-1:0]   m_axi_araddr,
    output wire [7:0]                   m_axi_arlen,
    output wire                         m_axi_arvalid,
    input  wire                         m_axi_arready,
    input  wire [AXI_DATA_WIDTH-1:0]   m_axi_rdata,
    input  wire                         m_axi_rlast,
    input  wire [1:0]                   m_axi_rresp,
    input  wire                         m_axi_rvalid,
    output wire                         m_axi_rready,
    
    // Interrupt
    output wire                         irq_frame_done,
    output wire                         irq_error
);

//==============================================================================
// Internal Signals
//==============================================================================

// Decoder State Machine
localparam STATE_IDLE         = 4'd0;
localparam STATE_PARSE_OBU    = 4'd1;
localparam STATE_PARSE_HEADER = 4'd2;
localparam STATE_DECODE_TILES = 4'd3;
localparam STATE_RECON        = 4'd4;
localparam STATE_FILTER       = 4'd5;
localparam STATE_OUTPUT       = 4'd6;
localparam STATE_ERROR        = 4'd7;

reg [3:0] decoder_state, decoder_state_next;

// OBU Parser Signals
wire        obu_valid;
wire [3:0]  obu_type;
wire [31:0] obu_size;
wire        obu_ready;
wire        obu_tready;

// Frame Header Signals
wire        frame_header_valid;
wire [1:0]  frame_type;  // 0:KEY, 1:INTER, 2:INTRA_ONLY, 3:S_FRAME
wire [15:0] frame_width;
wire [15:0] frame_height;
wire [7:0]  qindex;
wire        frame_header_ready;

// Tile Decoder Signals
wire        tile_valid;
wire [7:0]  tile_row;
wire [7:0]  tile_col;
wire        tile_ready;
wire        tile_done;

// Reconstruction Signals
wire        recon_valid;
wire        recon_ready;
wire        recon_done;

// Loop Filter Signals
wire        filter_valid;
wire        filter_ready;
wire        filter_done;

// Frame Buffer Signals - Internal Buses
wire [AXI_ADDR_WIDTH-1:0] tile_rd_addr;
wire                       tile_rd_en;
wire [AXI_ADDR_WIDTH-1:0] out_rd_addr;
wire                       out_rd_en;

wire [AXI_ADDR_WIDTH-1:0] fb_wr_addr;
wire [AXI_DATA_WIDTH-1:0] fb_wr_data;
wire                       fb_wr_en;
wire [AXI_ADDR_WIDTH-1:0] fb_rd_addr;
wire [AXI_DATA_WIDTH-1:0] fb_rd_data;
wire                       fb_rd_en;
wire [1:0]                fb_rd_sel_plane;  // 0: Y, 1: U, 2: V

// Status Registers
reg [31:0] frames_decoded;
reg [31:0] error_count;
reg        decoder_busy;

//==============================================================================
// Module Instantiations
//==============================================================================

// OBU Parser
av2_obu_parser #(
    .DATA_WIDTH(AXI_DATA_WIDTH)
) u_obu_parser (
    .clk            (clk),
    .rst_n          (rst_n),
    .s_axis_tdata   (s_axis_tdata),
    .s_axis_tvalid  (s_axis_tvalid),
    .s_axis_tready  (obu_tready),
    .s_axis_tlast   (s_axis_tlast),
    .obu_type       (obu_type),
    .obu_size       (obu_size),
    .obu_valid      (obu_valid),
    .obu_ready      (obu_ready)
);

// Frame Header Parser
av2_frame_header_parser u_header_parser (
    .clk                (clk),
    .rst_n              (rst_n),
    .obu_valid          (obu_valid && (obu_type == 4'd3)), // OBU_FRAME_HEADER
    .obu_data           (s_axis_tdata),
    .frame_type         (frame_type),
    .frame_width        (frame_width),
    .frame_height       (frame_height),
    .qindex             (qindex),
    .header_valid       (frame_header_valid),
    .header_ready       (decoder_state == STATE_PARSE_HEADER && frame_header_valid)
);

// Tile Decoder start pulse (only when entering STATE_DECODE_TILES)
reg tile_start_pulse_r;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n)
        tile_start_pulse_r <= 1'b0;
    else
        tile_start_pulse_r <= (decoder_state != STATE_DECODE_TILES) && (decoder_state_next == STATE_DECODE_TILES);
end
wire tile_start_pulse = tile_start_pulse_r;

av2_tile_decoder_complete #(
    .MAX_WIDTH  (MAX_FRAME_WIDTH),
    .MAX_HEIGHT (MAX_FRAME_HEIGHT),
    .PIXEL_WIDTH(PIXEL_WIDTH),
    .MAX_SB_SIZE(64)
) u_tile_decoder (
    .clk            (clk),
    .rst_n          (rst_n),
    .start          (tile_start_pulse),
    .frame_width    (frame_width),
    .frame_height   (frame_height),
    .qindex         (qindex),
    .frame_type     (frame_type),
    .tile_data      (s_axis_tdata),
    .tile_valid     (s_axis_tvalid),
    .tile_ready     (tile_ready),
    .ref_read_addr  (tile_rd_addr),
    .ref_pixel_data (fb_rd_data[9:0]),
    .ref_read_en    (tile_rd_en),
    .recon_data     (fb_wr_data),
    .recon_addr     (fb_wr_addr),
    .recon_wr_en    (fb_wr_en),
    .tile_done      (tile_done)
);

// Reconstruction and Loop Filter are integrated in av2_tile_decoder_complete module
assign recon_done = tile_done;
assign filter_done = tile_done;

// Frame Buffer Controller
av2_frame_buffer_ctrl #(
    .ADDR_WIDTH(AXI_ADDR_WIDTH),
    .DATA_WIDTH(AXI_DATA_WIDTH)
) u_fb_ctrl (
    .clk            (clk),
    .rst_n          (rst_n),
    // Internal Interface
    .wr_addr        (fb_wr_addr),
    .wr_data        (fb_wr_data),
    .wr_en          (fb_wr_en),
    .rd_addr        (fb_rd_addr),
    .rd_data        (fb_rd_data),
    .rd_en          (fb_rd_en),
    .rd_sel_plane   (fb_rd_sel_plane),
    // AXI4 Master Interface
    .m_axi_awaddr   (m_axi_awaddr),
    .m_axi_awlen    (m_axi_awlen),
    .m_axi_awvalid  (m_axi_awvalid),
    .m_axi_awready  (m_axi_awready),
    .m_axi_wdata    (m_axi_wdata),
    .m_axi_wlast    (m_axi_wlast),
    .m_axi_wvalid   (m_axi_wvalid),
    .m_axi_wready   (m_axi_wready),
    .m_axi_bresp    (m_axi_bresp),
    .m_axi_bvalid   (m_axi_bvalid),
    .m_axi_bready   (m_axi_bready),
    .m_axi_araddr   (m_axi_araddr),
    .m_axi_arlen    (m_axi_arlen),
    .m_axi_arvalid  (m_axi_arvalid),
    .m_axi_arready  (m_axi_arready),
    .m_axi_rdata    (m_axi_rdata),
    .m_axi_rlast    (m_axi_rlast),
    .m_axi_rresp    (m_axi_rresp),
    .m_axi_rvalid   (m_axi_rvalid),
    .m_axi_rready   (m_axi_rready)
);

// Output Controller start pulse (only when entering STATE_OUTPUT)
// Use a registered version to ensure single pulse
reg output_start_pulse_r;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n)
        output_start_pulse_r <= 1'b0;
    else
        output_start_pulse_r <= (decoder_state != STATE_OUTPUT) && (decoder_state_next == STATE_OUTPUT);
end
wire output_start_pulse = output_start_pulse_r;
wire output_done;

av2_output_ctrl #(
    .DATA_WIDTH(AXI_DATA_WIDTH)
) u_output_ctrl (
    .clk            (clk),
    .rst_n          (rst_n),
    .start          (output_start_pulse),
    .frame_width    (frame_width),
    .frame_height   (frame_height),
    .fb_rd_addr     (out_rd_addr),
    .fb_rd_sel_plane(fb_rd_sel_plane),
    .fb_rd_data     (fb_rd_data),
    .m_axis_tdata   (m_axis_tdata),
    .m_axis_tvalid  (m_axis_tvalid),
    .m_axis_tready  (m_axis_tready),
    .m_axis_tlast   (m_axis_tlast),
    .m_axis_tkeep   (m_axis_tkeep),
    .m_axis_tuser   (m_axis_tuser),
    .done           (output_done)
);

//==============================================================================
// Main Decoder State Machine
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n)
        decoder_state <= STATE_IDLE;
    else begin
        if (decoder_state != decoder_state_next) begin
            $display("[TIME %0t] Decoder State: %d -> %d", $time, decoder_state, decoder_state_next);
        end
        decoder_state <= decoder_state_next;
    end
end

always @(*) begin
    decoder_state_next = decoder_state;
    
    case (decoder_state)
        STATE_IDLE: begin
            if (s_axis_tvalid)
                decoder_state_next = STATE_PARSE_OBU;
        end
        
        STATE_PARSE_OBU: begin
            if (obu_valid) begin
                case (obu_type)
                    4'd1: decoder_state_next = STATE_IDLE;        // OBU_SEQUENCE_HEADER
                    4'd3: decoder_state_next = STATE_PARSE_HEADER; // OBU_FRAME_HEADER
                    4'd6: decoder_state_next = STATE_DECODE_TILES; // OBU_FRAME
                    default: decoder_state_next = STATE_IDLE;
                endcase
            end
        end
        
        STATE_PARSE_HEADER: begin
            if (frame_header_valid)
                decoder_state_next = STATE_DECODE_TILES;
        end
        
        STATE_DECODE_TILES: begin
            if (tile_done) begin
                decoder_state_next = STATE_OUTPUT;
            end
        end
        
        STATE_OUTPUT: begin
            // Wait for output controller to finish
            if (output_done)
                decoder_state_next = STATE_IDLE;
        end
        
        STATE_ERROR: begin
            decoder_state_next = STATE_IDLE;
        end
        
        default: decoder_state_next = STATE_IDLE;
    endcase
end

assign obu_ready = (decoder_state == STATE_PARSE_OBU);
assign frame_header_ready = (decoder_state == STATE_PARSE_HEADER);

//==============================================================================
// Control/Status Registers
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        frames_decoded <= 32'h0;
        error_count    <= 32'h0;
        decoder_busy   <= 1'b0;
    end else begin
        decoder_busy <= (decoder_state != STATE_IDLE);
        
        if (decoder_state == STATE_OUTPUT && decoder_state_next == STATE_IDLE)
            frames_decoded <= frames_decoded + 1'b1;
            
        if (decoder_state == STATE_ERROR)
            error_count <= error_count + 1'b1;
    end
end

// Register Read/Write
reg [31:0] reg_rdata_r;
reg        reg_ready_r;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        reg_rdata_r <= 32'h0;
        reg_ready_r <= 1'b0;
    end else begin
        reg_ready_r <= reg_rd_en || reg_wr_en;
        
        if (reg_rd_en) begin
            case (reg_addr)
                16'h0000: reg_rdata_r <= 32'h41563200;      // "AV2" signature
                16'h0004: reg_rdata_r <= {31'h0, decoder_busy};
                16'h0008: reg_rdata_r <= frames_decoded;
                16'h000C: reg_rdata_r <= error_count;
                16'h0010: reg_rdata_r <= {16'h0, frame_width};
                16'h0014: reg_rdata_r <= {16'h0, frame_height};
                16'h0018: reg_rdata_r <= {28'h0, decoder_state};
                default:  reg_rdata_r <= 32'hDEADBEEF;
            endcase
        end
    end
end

assign reg_ready = reg_rd_en || reg_wr_en;
assign reg_rdata = reg_rdata_r;

//==============================================================================
// Interrupt Generation
//==============================================================================

reg irq_frame_done_r;
reg irq_error_r;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        irq_frame_done_r <= 1'b0;
        irq_error_r      <= 1'b0;
    end else begin
        irq_frame_done_r <= (decoder_state == STATE_OUTPUT && decoder_state_next == STATE_IDLE);
        irq_error_r      <= (decoder_state == STATE_ERROR);
    end
end

assign irq_frame_done = irq_frame_done_r;
assign irq_error      = irq_error_r;

//==============================================================================
// Memory Interface Multiplexing
//==============================================================================

// Frame buffer read address multiplexing
assign fb_rd_addr = (decoder_state == STATE_OUTPUT) ? out_rd_addr : tile_rd_addr;

// Frame buffer read enable multiplexing
assign fb_rd_en   = (decoder_state == STATE_OUTPUT) ? 1'b1 : tile_rd_en;

// Input stream ready multiplexing
assign s_axis_tready = (decoder_state == STATE_DECODE_TILES) ? tile_ready : obu_tready;

endmodule
