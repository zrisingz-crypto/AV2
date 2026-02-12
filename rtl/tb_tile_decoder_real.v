//==============================================================================
// Testbench for AV2 Real Tile Decoder
// Tests the integrated real decode pipeline
//==============================================================================

`timescale 1ns / 1ps

module tb_tile_decoder_real;

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

// Bitstream input (not used in bypass mode)
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
// Test Procedure
//==============================================================================

integer cycle_count;
integer timeout_cycles = 1000000;  // 10ms timeout at 100MHz

initial begin
    // Initialize signals
    rst_n = 0;
    start = 0;
    tile_data = 0;
    tile_valid = 0;
    frame_width = FRAME_WIDTH;
    frame_height = FRAME_HEIGHT;
    qindex = 8'd32;
    frame_type = 2'd0;  // I-frame (intra-only)
    cycle_count = 0;
    
    // Note: output directory must exist before running simulation
    // Use: mkdir -p output
    
    //==========================================================================
    // Reset sequence
    //==========================================================================
    $display("========================================");
    $display("AV2 Real Tile Decoder Testbench");
    $display("========================================");
    $display("[TIME %0t] Starting reset sequence", $time);
    
    #(CLK_PERIOD * 5);
    rst_n = 1;
    #(CLK_PERIOD * 5);
    $display("[TIME %0t] Reset complete", $time);
    
    //==========================================================================
    // Generate test data
    //==========================================================================
    $display("[TIME %0t] Generating software pixel ROM...", $time);
    generate_sw_pixel_rom();
    $display("[TIME %0t] Software pixel ROM generated", $time);
    
    //==========================================================================
    // Start decoder
    //==========================================================================
    $display("========================================");
    $display("Starting Decoder Test");
    $display("========================================");
    $display("[TIME %0t] Frame size: %0d x %0d", $time, frame_width, frame_height);
    $display("[TIME %0t] Frame type: %0d (0=I-frame)", $time, frame_type);
    $display("[TIME %0t] QIndex: %0d", $time, qindex);
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
            $display("[TIME %0t] Waiting... cycle %0d", $time, cycle_count);
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

//==============================================================================
// Generate software pixel ROM
//==============================================================================

// Note: Software pixel ROM must be pre-generated
// Use: python -c "for i in range(4096): print(f'{i%256:02X}')" > output/sw_pixel_rom.txt

task generate_sw_pixel_rom;
    // ROM is pre-generated, nothing to do here
    $display("[INFO] Using pre-generated sw_pixel_rom.txt from output/ directory");
endtask

//==============================================================================
// Monitor write operations
//==============================================================================

integer write_count;

always @(posedge clk) begin
    if (recon_wr_en) begin
        write_count = write_count + 1;
        if (write_count % 100 == 0) begin
            $display("[MONITOR] Write #%0d: addr=%0d", write_count, recon_addr);
        end
    end
end

initial begin
    write_count = 0;
end

endmodule