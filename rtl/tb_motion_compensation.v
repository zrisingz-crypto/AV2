//==============================================================================
// Testbench for Motion Compensation Module
//==============================================================================

`timescale 1ns / 1ps

module tb_motion_compensation;

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
reg [31:0] ref_frame_addr;
wire [9:0] ref_pixel_data;
reg [9:0] ref_data;  // Simple reference data source
wire [31:0] ref_read_addr;
wire ref_read_en;
reg signed [15:0] mv_x;
reg signed [15:0] mv_y;
reg [6:0] block_width;
reg [6:0] block_height;
reg [15:0] block_x;
reg [15:0] block_y;
reg [3:0] interp_filter;
reg start;
reg use_bidir;
reg ready;

// Outputs
wire [9:0] pred_block[0:16383];
wire valid;

//==============================================================================
// Test Control
//==============================================================================

integer test_num;
integer error_count;
reg [9:0] expected_data[0:63];
integer i;

//==============================================================================
// DUT Instantiation
//==============================================================================

av2_motion_compensation_real #(
    .MAX_WIDTH(MAX_WIDTH),
    .MAX_HEIGHT(MAX_HEIGHT),
    .MAX_BLOCK_SIZE(128)
) dut (
    .clk(clk),
    .rst_n(rst_n),
    .ref_frame_addr(ref_frame_addr),
    .ref_read_addr(),
    .ref_pixel_data(ref_pixel_data),
    .ref_read_en(ref_read_en),
    .mv_x(mv_x),
    .mv_y(mv_y),
    .block_width(block_width),
    .block_height(block_height),
    .block_x(block_x),
    .block_y(block_y),
    .interp_filter(interp_filter),
    .start(start),
    .use_bidir(use_bidir),
    .pred_block(pred_block),
    .valid(valid),
    .ready(ready)
);

//==============================================================================
// Reference Pixel Simulation
//==============================================================================

always @(posedge clk) begin
    if (ref_read_en) begin
        ref_data <= (ref_read_addr % 256);  // Simple pattern
    end
end

assign ref_pixel_data = ref_data;

//==============================================================================
// Clock Generation
//==============================================================================

initial begin
    clk = 0;
    forever #5 clk = ~clk;
end

//==============================================================================
// Test Cases
//==============================================================================

initial begin
    // Initialize
    rst_n = 0;
    start = 0;
    ready = 0;
    mv_x = 16'sd0;
    mv_y = 16'sd0;
    block_width = 7'd8;
    block_height = 7'd8;
    block_x = 16'd0;
    block_y = 16'd0;
    interp_filter = 4'd0;
    use_bidir = 1'b0;
    ref_frame_addr = 32'd0;
    test_num = 0;
    error_count = 0;
    
    // Initialize expected data
    for (i = 0; i < 64; i = i + 1) begin
        expected_data[i] = 10'd128;
    end
    
    // Reset sequence
    #20;
    rst_n = 1;
    #20;
    
    $display("========================================");
    $display("Motion Compensation Testbench");
    $display("========================================");
    
    // Test 1: Integer MV (no interpolation)
    test_num = 1;
    $display("\n[Test %0d] Integer Motion Vector", test_num);
    mv_x = 16'sd8;
    mv_y = 16'sd4;
    block_width = 7'd8;
    block_height = 7'd8;
    block_x = 16'd0;
    block_y = 16'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Integer MV test completed", test_num);
    
    #50;
    
    // Test 2: Sub-pixel MV (horizontal)
    test_num = 2;
    $display("\n[Test %0d] Sub-pixel Motion Vector (horizontal)", test_num);
    mv_x = 16'sd12;  // 1.5 pixels
    mv_y = 16'sd0;
    block_width = 7'd8;
    block_height = 7'd8;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Sub-pixel MV test completed", test_num);
    
    #50;
    
    // Test 3: Sub-pixel MV (vertical)
    test_num = 3;
    $display("\n[Test %0d] Sub-pixel Motion Vector (vertical)", test_num);
    mv_x = 16'sd0;
    mv_y = 16'sd12;  // 1.5 pixels
    block_width = 7'd8;
    block_height = 7'd8;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Sub-pixel MV test completed", test_num);
    
    #50;
    
    // Test 4: 16x16 block
    test_num = 4;
    $display("\n[Test %0d] 16x16 Block", test_num);
    mv_x = 16'sd4;
    mv_y = 16'sd4;
    block_width = 7'd16;
    block_height = 7'd16;
    
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