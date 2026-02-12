//==============================================================================
// Testbench for AV2 Real IVF Decoding
// Tests the integrated real decode pipeline with actual IVF bitstream
//==============================================================================

`timescale 1ns / 1ps

module tb_ivf_decode_real;

//==============================================================================
// Parameters
//==============================================================================

parameter CLK_PERIOD = 10;  // 100 MHz clock

parameter FRAME_WIDTH  = 64;
parameter FRAME_HEIGHT = 64;
parameter MAX_SB_SIZE  = 16;

//==============================================================================
// Signals
//==============================================================================

// Clock and reset
reg clk;
reg rst_n;

// Control
reg start;

// Frame parameters
reg [15:0] frame_width;
reg [15:0] frame_height;
reg [7:0]  qindex;
reg [1:0]  frame_type;

// Bitstream input
reg [127:0] tile_data;
reg        tile_valid;
wire       tile_ready;

// Reference frame interface (not used in intra-only mode)
wire [31:0] ref_read_addr;
reg  [9:0]  ref_pixel_data;
wire        ref_read_en;

// Reconstructed frame output
wire [127:0] recon_data;
wire [31:0]  recon_addr;
wire        recon_wr_en;
wire        tile_done;

//==============================================================================
// DUT Instantiation
//==============================================================================

av2_tile_decoder_real #(
    .MAX_WIDTH(FRAME_WIDTH),
    .MAX_HEIGHT(FRAME_HEIGHT),
    .PIXEL_WIDTH(10),
    .MAX_SB_SIZE(MAX_SB_SIZE)
) dut (
    .clk              (clk),
    .rst_n            (rst_n),
    .start            (start),
    
    // Frame parameters
    .frame_width      (frame_width),
    .frame_height     (frame_height),
    .qindex           (qindex),
    .frame_type       (frame_type),
    
    // Bitstream input
    .tile_data        (tile_data),
    .tile_valid       (tile_valid),
    .tile_ready       (tile_ready),
    
    // Reference frame interface
    .ref_read_addr    (ref_read_addr),
    .ref_pixel_data   (ref_pixel_data),
    .ref_read_en      (ref_read_en),
    
    // Reconstructed frame output
    .recon_data       (recon_data),
    .recon_addr       (recon_addr),
    .recon_wr_en      (recon_wr_en),
    .tile_done        (tile_done)
);

//==============================================================================
// Clock Generation
//==============================================================================

initial begin
    clk = 0;
    forever #(CLK_PERIOD/2) clk = ~clk;
end

//==============================================================================
// Reference pixel data (dummy, not used in intra-only mode)
//==============================================================================

always @(*) begin
    ref_pixel_data = 10'd128;
end

//==============================================================================
// IVF Bitstream Reader
//==============================================================================

// Instantiate IVF bitstream data module
wire [7:0] frame_0_data_wire [0:3391];
wire [7:0] frame_1_data_wire [0:46];

ivf_bitstream_data #(
    .FRAME_0_SIZE (3392),
    .FRAME_1_SIZE (47),
    .IVF_NUM_FRAMES (2)
) u_bitstream_data (
    .frame_0_data (frame_0_data_wire),
    .frame_1_data (frame_1_data_wire)
);

reg [31:0] bitstream_ptr;      // Current position in bitstream (bytes)
reg [31:0] bitstream_bytes;    // Total bytes in current frame
reg [7:0]  current_frame_data[0:4095]; // Buffer for current frame
reg        frame_data_ready;
reg        reading_frame;

// Frame counter
reg [31:0] frame_num;
reg [31:0] total_frames;

initial begin
    frame_num = 0;
    total_frames = 2;
    bitstream_ptr = 0;
    bitstream_bytes = 0;
    frame_data_ready = 0;
    reading_frame = 0;
    
    // Load first frame data
    load_frame_data(0);
end

task load_frame_data;
    input [31:0] frame_idx;
    integer i;
    begin
        if (frame_idx == 0) begin
            // Load frame 0
            bitstream_bytes = 3392;
            for (i = 0; i < 3392; i = i + 1) begin
                current_frame_data[i] = frame_0_data_wire[i];
            end
        end else if (frame_idx == 1) begin
            // Load frame 1
            bitstream_bytes = 47;
            for (i = 0; i < 47; i = i + 1) begin
                current_frame_data[i] = frame_1_data_wire[i];
            end
        end
        
        bitstream_ptr = 0;
        frame_data_ready = 1;
        $display("[TIME %0t] Loaded frame %0d: %0d bytes", $time, frame_idx, bitstream_bytes);
    end
endtask

//==============================================================================
// Bitstream output to DUT
//==============================================================================

reg [127:0] tile_data_reg;
reg        tile_valid_reg;
integer    words_sent;

always @(posedge clk) begin
    if (rst_n && tile_ready && frame_data_ready) begin
        // Pack 16 bytes into 128-bit word
        tile_data_reg = {current_frame_data[bitstream_ptr+15],
                       current_frame_data[bitstream_ptr+14],
                       current_frame_data[bitstream_ptr+13],
                       current_frame_data[bitstream_ptr+12],
                       current_frame_data[bitstream_ptr+11],
                       current_frame_data[bitstream_ptr+10],
                       current_frame_data[bitstream_ptr+9],
                       current_frame_data[bitstream_ptr+8],
                       current_frame_data[bitstream_ptr+7],
                       current_frame_data[bitstream_ptr+6],
                       current_frame_data[bitstream_ptr+5],
                       current_frame_data[bitstream_ptr+4],
                       current_frame_data[bitstream_ptr+3],
                       current_frame_data[bitstream_ptr+2],
                       current_frame_data[bitstream_ptr+1],
                       current_frame_data[bitstream_ptr+0]};
        
        tile_valid_reg = 1'b1;
        
        if (bitstream_ptr + 16 < bitstream_bytes) begin
            bitstream_ptr <= bitstream_ptr + 16;
        end else begin
            // End of frame data
            bitstream_ptr <= 32'd0;
            tile_valid_reg = 1'b0;
        end
        
        words_sent = words_sent + 1;
        if (words_sent < 10) begin
            $display("[TIME %0t] Sent bitstream word %0d: %032x", $time, words_sent, tile_data_reg);
        end
    end else begin
        tile_valid_reg = 1'b0;
    end
end

assign tile_data = tile_data_reg;
assign tile_valid = tile_valid_reg;

//==============================================================================
// Output File Writing
//==============================================================================

integer output_file;
integer write_count;

initial begin
    output_file = $fopen("output/rtl_decoded_frame0.yuv", "w");
    $display("[INIT] fopen returned handle=%0d", output_file);
    if (output_file == 0) begin
        $display("ERROR: Cannot open output file");
        $finish;
    end
    write_count = 0;
    $display("[INIT] Output file opened successfully, file handle=%0d", output_file);
end

always @(posedge clk) begin
    if (recon_wr_en) begin
        // Write 128-bit word (16 bytes) to file
        $fwrite(output_file, "%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c",
                recon_data[7:0], recon_data[15:8], recon_data[23:16], recon_data[31:24],
                recon_data[39:32], recon_data[47:40], recon_data[55:48], recon_data[63:56],
                recon_data[71:64], recon_data[79:72], recon_data[87:80], recon_data[95:88],
                recon_data[103:96], recon_data[111:104], recon_data[119:112], recon_data[127:120]);
        write_count = write_count + 1;
        
        // Explicit flush every write to ensure data is written
        $fflush(output_file);
        
        if (write_count == 1) begin
            $display("[FILE] First write: count=%0d, recon_data=%016x", write_count, recon_data);
        end
        
        if (write_count % 64 == 0) begin
            $display("[MONITOR] Write #%0d: addr=%0d", $time, recon_addr);
            $display("[TIME %0t] Written %0d bytes (%0d words)", $time, write_count * 16, write_count);
        end
    end else if (write_count > 0 && write_count % 100 == 0) begin
        // Periodic debug to show we're monitoring
        $display("[MONITOR] No write, cycle count: %0d", $time);
    end
end

// Debug: Show when recon_wr_en is asserted
always @(posedge clk) begin
    if (recon_wr_en) begin
        $display("[DEBUG] recon_wr_en=1 at time %0t, recon_data=%016x", $time, recon_data);
    end
end

final begin
    $fclose(output_file);
    $display("[INFO] Total bytes written: %0d", write_count * 16);
end

//==============================================================================
// Test Procedure
//==============================================================================

integer cycle_count;
integer timeout_cycles = 1000000;  // 10ms timeout at 100MHz

initial begin
    // Initialize signals
    rst_n = 0;
    start = 0;
    frame_width = FRAME_WIDTH;
    frame_height = FRAME_HEIGHT;
    qindex = 8'd32;
    frame_type = 2'd0;  // I-frame (intra-only)
    cycle_count = 0;
    words_sent = 0;
    
    //==========================================================================
    // Reset sequence
    //==========================================================================
    $display("========================================");
    $display("AV2 Real IVF Decoder Testbench");
    $display("========================================");
    $display("[TIME %0t] Starting reset sequence", $time);
    
    #(CLK_PERIOD * 5);
    rst_n = 1;
    #(CLK_PERIOD * 5);
    $display("[TIME %0t] Reset complete", $time);
    
    //==========================================================================
    // Start decoder
    //==========================================================================
    $display("========================================");
    $display("Starting IVF Decoder Test");
    $display("========================================");
    $display("[TIME %0t] Frame size: %0d x %0d", $time, frame_width, frame_height);
    $display("[TIME %0t] Frame type: %0d (0=I-frame)", $time, frame_type);
    $display("[TIME %0t] QIndex: %0d", $time, qindex);
    $display("[TIME %0t] Bitstream size: %0d bytes", $time, bitstream_bytes);
    $display("[TIME %0t] Starting decoder...", $time);
    
    #(CLK_PERIOD);
    start = 1;
    #(CLK_PERIOD);
    start = 0;
    $display("[TIME %0t] Start signal asserted", $time);
    
    //==========================================================================
    // Wait for completion or timeout
    //==========================================================================
    $display("[TIME %0t] Waiting for tile_done signal...", $time);
    
    while (!tile_done && cycle_count < timeout_cycles) begin
        #(CLK_PERIOD);
        cycle_count = cycle_count + 1;
        
        // Print progress every 1000 cycles
        if (cycle_count % 1000 == 0) begin
            $display("[TIME %0t] Waiting... cycle %0d, words_sent=%0d", $time, cycle_count, words_sent);
        end
    end
    
    //==========================================================================
    // Check results
    //==========================================================================
    if (tile_done) begin
        $display("========================================");
        $display("TEST PASSED!");
        $display("========================================");
        $display("[TIME %0t] tile_done received at cycle %0d", $time, cycle_count);
        $display("[TIME %0t] Total cycles: %0d", $time, cycle_count);
        $display("[TIME %0t] Simulation time: %0d ns", $time, $time);
        $display("[TIME %0t] Output written to output/rtl_decoded_frame0.yuv", $time);
        #(CLK_PERIOD * 100);
        $finish;
    end else begin
        $display("========================================");
        $display("TEST FAILED - TIMEOUT!");
        $display("========================================");
        $display("[TIME %0t] Timeout after %0d cycles", $time, cycle_count);
        $display("[TIME %0t] Simulation time: %0d ns", $time, $time);
        $finish;
    end
end

endmodule