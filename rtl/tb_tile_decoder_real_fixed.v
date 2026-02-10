//==============================================================================
// Testbench for AV2 Real Tile Decoder (FIXED VERSION)
// Tests the integrated real decode modules
//==============================================================================

`timescale 1ns / 1ps

module tb_tile_decoder_real_fixed;

//==============================================================================
// Parameters
//==============================================================================

parameter CLK_PERIOD = 10;  // 100MHz
parameter MAX_WIDTH = 64;
parameter MAX_HEIGHT = 64;
parameter PIXEL_WIDTH = 10;
parameter TIMEOUT_CYCLES = 5000;  // Increased timeout for real modules

//==============================================================================
// Signals
//==============================================================================

reg clk;
reg rst_n;
reg start;

// Frame parameters
reg [15:0] frame_width;
reg [15:0] frame_height;
reg [7:0]  qindex;
reg [1:0]  frame_type;

// Bitstream input
reg [127:0] tile_data;
reg         tile_valid;
wire        tile_ready;

// Reference frame interface (for inter prediction)
wire [31:0] ref_read_addr;
reg  [9:0]  ref_pixel_data;
wire        ref_read_en;

// Reconstructed frame output
wire [127:0] recon_data;
wire [31:0]  recon_addr;
wire        recon_wr_en;
wire        tile_done;

// Internal signals for observation
integer cycle_count;
integer decode_cycle;

// Output file
integer rtl_output_file;

//==============================================================================
// DUT Instance
//==============================================================================

av2_tile_decoder_real_fixed #(
    .MAX_WIDTH(MAX_WIDTH),
    .MAX_HEIGHT(MAX_HEIGHT),
    .PIXEL_WIDTH(PIXEL_WIDTH)
) dut (
    .clk              (clk),
    .rst_n            (rst_n),
    .start            (start),
    .frame_width      (frame_width),
    .frame_height     (frame_height),
    .qindex           (qindex),
    .frame_type       (frame_type),
    .tile_data        (tile_data),
    .tile_valid       (tile_valid),
    .tile_ready       (tile_ready),
    .ref_read_addr    (ref_read_addr),
    .ref_pixel_data   (ref_pixel_data),
    .ref_read_en      (ref_read_en),
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
// Reset Generation
//==============================================================================

initial begin
    rst_n = 0;
    #(CLK_PERIOD * 5);
    rst_n = 1;
end

//==============================================================================
// Test Process
//==============================================================================

integer i;

initial begin
    // Open output file
    rtl_output_file = $fopen("rtl_real_output.txt", "w");
    if (rtl_output_file == 0) begin
        $display("ERROR: Cannot open output file!");
        $finish;
    end
    
    // Initialize signals
    start = 0;
    tile_data = 128'h0;
    tile_valid = 0;
    ref_pixel_data = 10'd0;
    frame_width = 16'd64;
    frame_height = 16'd64;
    qindex = 8'd128;
    frame_type = 2'd0;  // Key frame (I-frame)
    
    cycle_count = 0;
    decode_cycle = 0;
    
    // Wait for reset
    @(posedge rst_n);
    @(posedge clk);
    
    $display("\n========================================");
    $display("Test 1: 64x64 I-frame decoding");
    $display("========================================\n");
    
    // Start decoding
    $display("[TIME %0t] Starting decode...", $time);
    start <= 1;
    @(posedge clk);
    start <= 0;
    
    // Provide bitstream data (simplified: use varying patterns)
    for (i = 0; i < 100; i = i + 1) begin
        @(posedge clk);
        
        if (tile_ready) begin
            tile_valid = 1;
            // Generate test bitstream data - use different patterns
            tile_data = {i[15:0], i[15:0], i[15:0], i[15:0], 
                         i[15:0], i[15:0], i[15:0], i[15:0]};
        end else begin
            tile_valid = 0;
        end
    end
    
    tile_valid = 0;
    
    // Wait for decode to complete
    $display("[TIME %0t] Waiting for tile_done...", $time);
    @(posedge tile_done);
    $display("[TIME %0t] Decode complete! Total cycles: %0d", $time, cycle_count);
    
    // Wait a bit more for output writes to complete
    #(CLK_PERIOD * 100);
    
    // Close output file
    $fclose(rtl_output_file);
    
    $display("\n========================================");
    $display("Test Complete!");
    $display("========================================");
    $display("Total cycles: %0d", cycle_count);
    $display("Output file: rtl_real_output.txt");
    $display("========================================\n");
    
    #(CLK_PERIOD * 10);
    $finish;
end

//==============================================================================
// Monitor
//==============================================================================

// Monitor cycle count
always @(posedge clk) begin
    if (!rst_n)
        cycle_count = 0;
    else
        cycle_count = cycle_count + 1;
end

// Monitor reconstruction output
always @(posedge clk) begin
    if (recon_wr_en) begin
        integer j;
        reg [7:0] pixel;
        
        for (j = 0; j < 16; j = j + 1) begin
            case (j)
                4'd0:  pixel = recon_data[7:0];
                4'd1:  pixel = recon_data[15:8];
                4'd2:  pixel = recon_data[23:16];
                4'd3:  pixel = recon_data[31:24];
                4'd4:  pixel = recon_data[39:32];
                4'd5:  pixel = recon_data[47:40];
                4'd6:  pixel = recon_data[55:48];
                4'd7:  pixel = recon_data[63:56];
                4'd8:  pixel = recon_data[71:64];
                4'd9:  pixel = recon_data[79:72];
                4'd10: pixel = recon_data[87:80];
                4'd11: pixel = recon_data[95:88];
                4'd12: pixel = recon_data[103:96];
                4'd13: pixel = recon_data[111:104];
                4'd14: pixel = recon_data[119:112];
                4'd15: pixel = recon_data[127:120];
            endcase
            
            // Write to file
            $fwrite(rtl_output_file, "pixel[%0d][%0d] = %0d\n", 
                     recon_addr, j, pixel);
        end
        
        // Print progress every 256 bytes
        if (recon_addr % 256 == 0)
            $display("[TIME %0t] Write: addr=%0d, data=%0h", 
                     $time, recon_addr, recon_data);
    end
end

// Monitor reference frame reads
always @(posedge clk) begin
    if (ref_read_en) begin
        // Provide dummy reference pixel data
        ref_pixel_data = 10'd128;
    end
end

//==============================================================================
// Monitor for completion
//==============================================================================

// Monitor that tile_done is asserted
reg [31:0] start_cycle;
reg        start_seen;

always @(posedge clk) begin
    if (!rst_n) begin
        start_seen <= 1'b0;
    end else if (start) begin
        start_cycle <= cycle_count;
        start_seen <= 1'b1;
        $display("[TIME %0t] Start detected, cycle=%0d", $time, cycle_count);
    end
end

// Check timeout
always @(posedge clk) begin
    if (rst_n && start_seen && !tile_done) begin
        if (cycle_count - start_cycle > TIMEOUT_CYCLES) begin
            $display("[ERROR] tile_done not asserted within %0d cycles! Current cycle: %0d", 
                     TIMEOUT_CYCLES, cycle_count);
            $display("[ERROR] Possible deadlock - check state machine");
            $fclose(rtl_output_file);
            $finish;
        end
    end
end

// Print statistics
final begin
    $display("\n========================================");
    $display("Statistics");
    $display("========================================");
    $display("Total cycles: %0d", cycle_count);
    $display("========================================\n");
end

endmodule
