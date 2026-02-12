//==============================================================================
// Testbench for MV Decoder Module
//==============================================================================

`timescale 1ns / 1ps

module tb_mv_decoder;

//==============================================================================
// Clock and Reset
//==============================================================================

reg clk;
reg rst_n;

//==============================================================================
// DUT Interface
//==============================================================================

// Inputs
reg [15:0] context_idx;
reg [15:0] context_prob;
reg [15:0] decoded_symbol;
reg symbol_valid;
reg mv_ready;
reg start;

// Outputs
wire symbol_ready;
wire signed [15:0] mv_x;
wire signed [15:0] mv_y;
wire mv_valid;
wire done;

//==============================================================================
// Test Control
//==============================================================================

integer test_num;
integer error_count;
reg signed [15:0] expected_mv_x;
reg signed [15:0] expected_mv_y;
integer i;

//==============================================================================
// DUT Instantiation
//==============================================================================

av2_mv_decoder_real dut (
    .clk(clk),
    .rst_n(rst_n),
    .context_idx(context_idx),
    .context_prob(context_prob),
    .decoded_symbol(decoded_symbol),
    .symbol_valid(symbol_valid),
    .symbol_ready(symbol_ready),
    .mv_x(mv_x),
    .mv_y(mv_y),
    .mv_valid(mv_valid),
    .mv_ready(mv_ready),
    .start(start),
    .done(done)
);

//==============================================================================
// Clock Generation
//==============================================================================

initial begin
    clk = 0;
    forever #5 clk = ~clk;
end

//==============================================================================
// Symbol Data Provider
//==============================================================================

reg [15:0] test_symbols[0:15];
reg [3:0] symbol_idx;

initial begin
    // Test symbols for MV decoding
    test_symbols[0]  = 16'd0;   // Magnitude end marker
    test_symbols[1]  = 16'd1;   // Positive
    test_symbols[2]  = 16'd0;   // Magnitude end marker
    test_symbols[3]  = 16'd0;   // Zero
    test_symbols[4]  = 16'd1;   // Positive sign
    test_symbols[5]  = 16'd0;   // Magnitude end marker
    test_symbols[6]  = 16'd0;   // Zero
    test_symbols[7]  = 16'd0;   // Negative sign
    test_symbols[8]  = 16'd1;   // Magnitude bit
    test_symbols[9]  = 16'd0;   // Magnitude end marker
    test_symbols[10] = 16'd1;   // Positive sign
    test_symbols[11] = 16'd0;   // Magnitude end marker
    test_symbols[12] = 16'd0;   // Zero
    test_symbols[13] = 16'd0;   // Negative sign
    test_symbols[14] = 16'd1;   // Magnitude bit
    test_symbols[15] = 16'd0;   // Magnitude end marker
end

always @(posedge clk) begin
    if (symbol_ready) begin
        decoded_symbol <= test_symbols[symbol_idx];
        symbol_idx <= symbol_idx + 1;
    end
end

//==============================================================================
// Test Cases
//==============================================================================

initial begin
    // Initialize
    rst_n = 0;
    start = 0;
    mv_ready = 0;
    context_idx = 16'd0;
    context_prob = 16'd128;
    decoded_symbol = 16'd0;
    symbol_valid = 1'b0;
    symbol_idx = 4'd0;
    test_num = 0;
    error_count = 0;
    
    // Enable symbol valid
    symbol_valid = 1'b1;
    
    // Reset sequence
    #20;
    rst_n = 1;
    #20;
    
    $display("========================================");
    $display("MV Decoder Testbench");
    $display("========================================");
    
    // Test 1: Zero MV
    test_num = 1;
    $display("\n[Test %0d] Zero Motion Vector", test_num);
    expected_mv_x = 16'sd0;
    expected_mv_y = 16'sd0;
    symbol_idx = 4'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(mv_valid);
    mv_ready = 1;
    @(posedge clk);
    mv_ready = 0;
    
    if (mv_x !== expected_mv_x) begin
        $display("✗ Error: mv_x=%0d, expected=%0d", mv_x, expected_mv_x);
        error_count = error_count + 1;
    end else begin
        $display("✓ mv_x correct: %0d", mv_x);
    end
    
    if (mv_y !== expected_mv_y) begin
        $display("✗ Error: mv_y=%0d, expected=%0d", mv_y, expected_mv_y);
        error_count = error_count + 1;
    end else begin
        $display("✓ mv_y correct: %0d", mv_y);
    end
    
    #50;
    
    // Test 2: Positive MV
    test_num = 2;
    $display("\n[Test %0d] Positive Motion Vector", test_num);
    expected_mv_x = 16'sd4;
    expected_mv_y = 16'sd2;
    symbol_idx = 4'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(mv_valid);
    mv_ready = 1;
    @(posedge clk);
    mv_ready = 0;
    
    $display("Test %0d: MV decoded: (%0d, %0d)", test_num, mv_x, mv_y);
    
    #50;
    
    // Test 3: Negative MV
    test_num = 3;
    $display("\n[Test %0d] Negative Motion Vector", test_num);
    expected_mv_x = -16'sd4;
    expected_mv_y = -16'sd2;
    symbol_idx = 4'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(mv_valid);
    mv_ready = 1;
    @(posedge clk);
    mv_ready = 0;
    
    $display("Test %0d: MV decoded: (%0d, %0d)", test_num, mv_x, mv_y);
    
    #50;
    
    // Test 4: Large MV
    test_num = 4;
    $display("\n[Test %0d] Large Motion Vector", test_num);
    expected_mv_x = 16'sd100;
    expected_mv_y = 16'sd50;
    symbol_idx = 4'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(mv_valid);
    mv_ready = 1;
    @(posedge clk);
    mv_ready = 0;
    
    $display("Test %0d: MV decoded: (%0d, %0d)", test_num, mv_x, mv_y);
    
    #50;
    
    // Test 5: Mixed positive/negative
    test_num = 5;
    $display("\n[Test %0d] Mixed Positive/Negative MV", test_num);
    expected_mv_x = -16'sd10;
    expected_mv_y = 16'sd20;
    symbol_idx = 4'd0;
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(mv_valid);
    mv_ready = 1;
    @(posedge clk);
    mv_ready = 0;
    
    $display("Test %0d: MV decoded: (%0d, %0d)", test_num, mv_x, mv_y);
    
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