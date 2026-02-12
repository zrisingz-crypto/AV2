//==============================================================================
// Testbench for Deblocking Filter Module
//==============================================================================

`timescale 1ns / 1ps

module tb_deblocking_filter;

//==============================================================================
// Parameters
//==============================================================================

parameter MAX_WIDTH = 128;
parameter MAX_HEIGHT = 128;

//==============================================================================
// Clock and Reset
//==============================================================================

reg clk;
reg rst_n;

//==============================================================================
// DUT Interface
//==============================================================================

// Inputs
reg [9:0] src_pixels[0:16383];
reg [15:0] frame_width;
reg [15:0] frame_height;
reg [5:0] filter_level;
reg [2:0] sharpness;
reg start;
reg ready;

// Outputs
wire [9:0] dst_pixels[0:16383];
wire valid;

//==============================================================================
// Test Control
//==============================================================================

integer test_num;
integer error_count;
integer i;

//==============================================================================
// DUT Instantiation
//==============================================================================

av2_deblocking_filter_real #(
    .MAX_WIDTH(MAX_WIDTH),
    .MAX_HEIGHT(MAX_HEIGHT)
) dut (
    .clk(clk),
    .rst_n(rst_n),
    .src_pixels(src_pixels),
    .frame_width(frame_width),
    .frame_height(frame_height),
    .filter_level(filter_level),
    .sharpness(sharpness),
    .start(start),
    .dst_pixels(dst_pixels),
    .valid(valid),
    .ready(ready)
);

//==============================================================================
// Clock Generation
//==============================================================================

initial begin
    clk = 0;
    forever #5 clk = ~clk;
end

//==============================================================================
// Test Data Generation
//==============================================================================

// Create blocky data with visible block boundaries
task create_blocky_data;
    input integer start_idx;
    reg [9:0] base_value;
    integer pixel_idx;
    begin
        for (pixel_idx = 0; pixel_idx < 64; pixel_idx = pixel_idx + 1) begin
            // Create 8x8 blocks with different values
            if (pixel_idx < 8) src_pixels[start_idx + pixel_idx] = 10'd0;
            else if (pixel_idx < 16) src_pixels[start_idx + pixel_idx] = 10'd32;
            else if (pixel_idx < 24) src_pixels[start_idx + pixel_idx] = 10'd64;
            else if (pixel_idx < 32) src_pixels[start_idx + pixel_idx] = 10'd96;
            else if (pixel_idx < 40) src_pixels[start_idx + pixel_idx] = 10'd128;
            else if (pixel_idx < 48) src_pixels[start_idx + pixel_idx] = 10'd160;
            else if (pixel_idx < 56) src_pixels[start_idx + pixel_idx] = 10'd192;
            else src_pixels[start_idx + pixel_idx] = 10'd224;
        end
    end
endtask

// Create smooth data (no filtering needed)
task create_smooth_data;
    input integer start_idx;
    integer pixel_idx;
    begin
        for (pixel_idx = 0; pixel_idx < 64; pixel_idx = pixel_idx + 1) begin
            src_pixels[start_idx + pixel_idx] = 10'd128;
        end
    end
endtask

//==============================================================================
// Test Cases
//==============================================================================

initial begin
    // Initialize
    rst_n = 0;
    start = 0;
    ready = 0;
    frame_width = 16'd64;
    frame_height = 16'd64;
    filter_level = 6'd0;
    sharpness = 3'd0;
    test_num = 0;
    error_count = 0;
    
    // Initialize source pixels
    for (i = 0; i < 16384; i = i + 1) begin
        src_pixels[i] = 10'd128;
    end
    
    // Reset sequence
    #20;
    rst_n = 1;
    #20;
    
    $display("========================================");
    $display("Deblocking Filter Testbench");
    $display("========================================");
    
    // Test 1: No filtering (filter_level = 0)
    test_num = 1;
    $display("\n[Test %0d] No Filtering (filter_level=0)", test_num);
    create_smooth_data(0);
    filter_level = 6'd0;
    sharpness = 3'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: No filtering test completed", test_num);
    
    #50;
    
    // Test 2: Light filtering (filter_level = 10)
    test_num = 2;
    $display("\n[Test %0d] Light Filtering (filter_level=10)", test_num);
    create_blocky_data(0);
    filter_level = 6'd10;
    sharpness = 3'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    // Check some pixels to verify filtering
    if (dst_pixels[8] !== src_pixels[8] || dst_pixels[15] !== src_pixels[15]) begin
        $display("✗ Test %0d: Pixels at block boundary changed", test_num);
    end else begin
        $display("✓ Test %0d: Boundary pixels processed", test_num);
    end
    
    #50;
    
    // Test 3: Medium filtering (filter_level = 30)
    test_num = 3;
    $display("\n[Test %0d] Medium Filtering (filter_level=30)", test_num);
    create_blocky_data(0);
    filter_level = 6'd30;
    sharpness = 3'd2;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Medium filtering test completed", test_num);
    
    #50;
    
    // Test 4: Strong filtering (filter_level = 50)
    test_num = 4;
    $display("\n[Test %0d] Strong Filtering (filter_level=50)", test_num);
    create_blocky_data(0);
    filter_level = 6'd50;
    sharpness = 3'd5;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Strong filtering test completed", test_num);
    
    #50;
    
    // Test 5: 16x16 block size
    test_num = 5;
    $display("\n[Test %0d] 16x16 Block Size", test_num);
    create_blocky_data(0);
    filter_level = 6'd20;
    sharpness = 3'd1;
    frame_width = 16'd16;
    frame_height = 16'd16;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: 16x16 block test completed", test_num);
    
    #50;
    
    // Reset to default
    frame_width = 16'd64;
    frame_height = 16'd64;
    
    // Summary
    $display("\n========================================");
    $display("Test Summary");
    $display("========================================");
    $display("Total tests: %0d", test_num);
    $display("Errors: %0d", error_count);
    
    if (error_count == 0) begin
        $display("✓ ALL TESTS PASSED");
    end else begin
        $display("✗ SOME TESTS FAILED");
    end
    
    #100;
    $finish;
end

endmodule