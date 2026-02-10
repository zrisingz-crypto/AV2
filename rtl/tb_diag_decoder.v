//==============================================================================
// Testbench for Diagnostic Decoder
//==============================================================================

`timescale 1ns / 1ps

module tb_diag_decoder;

parameter CLK_PERIOD = 10;

reg clk;
reg rst_n;
reg start;
reg [15:0] frame_width;
reg [15:0] frame_height;
wire tile_done;

// DUT
diag_tile_decoder dut (
    .clk(clk),
    .rst_n(rst_n),
    .start(start),
    .frame_width(frame_width),
    .frame_height(frame_height),
    .tile_done(tile_done)
);

// Clock
initial begin
    clk = 0;
    forever #(CLK_PERIOD/2) clk = ~clk;
end

// Reset
initial begin
    rst_n = 0;
    #(CLK_PERIOD * 5);
    rst_n = 1;
end

// Test
initial begin
    start = 0;
    frame_width = 64;
    frame_height = 64;
    
    @(posedge rst_n);
    @(posedge clk);
    @(posedge clk);
    
    $display("\n========================================");
    $display("Diagnostic Decoder Test");
    $display("========================================\n");
    
    @(posedge clk);
    $display("[TIME %0t] Starting decode...", $time);
    start <= 1;  // Use non-blocking assignment
    @(posedge clk);
    start <= 0;
    
    // Wait for completion
    wait(tile_done);
    @(posedge clk);
    $display("[TIME %0t] Decode complete!", $time);
    
    @(posedge clk);
    @(posedge clk);
    
    $display("\n========================================");
    $display("Test PASSED!");
    $display("========================================\n");
    
    $finish;
end

// Timeout check
initial begin
    #(CLK_PERIOD * 1000);
    $display("[ERROR] Timeout!");
    $finish;
end

endmodule
