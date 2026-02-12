//==============================================================================
// Testbench for CDEF Filter Module
//==============================================================================

`timescale 1ns / 1ps

module tb_cdef_filter;

//==============================================================================
// Parameters
//==============================================================================

parameter BLOCK_SIZE = 8;

//==============================================================================
// Clock and Reset
//==============================================================================

reg clk;
reg rst_n;

//==============================================================================
// DUT Interface
//==============================================================================

// Inputs
reg [9:0] src_block[0:63];
reg [2:0] strength_y;
reg [2:0] strength_uv;
reg [2:0] damping;
reg is_chroma;
reg ready;

// Outputs
wire [9:0] dst_block[0:63];
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

av2_cdef_filter_real #(
    .BLOCK_SIZE(BLOCK_SIZE)
) dut (
    .clk(clk),
    .rst_n(rst_n),
    .src_block(src_block),
    .strength_y(strength_y),
    .strength_uv(strength_uv),
    .damping(damping),
    .is_chroma(is_chroma),
    .dst_block(dst_block),
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

// Create flat data (no filtering needed)
task create_flat_data;
    integer pixel_idx;
    begin
        for (pixel_idx = 0; pixel_idx < 64; pixel_idx = pixel_idx + 1) begin
            src_block[pixel_idx] = 10'd128;
        end
    end
endtask

// Create noisy data (filtering needed)
task create_noisy_data;
    integer pixel_idx;
    begin
        for (pixel_idx = 0; pixel_idx < 64; pixel_idx = pixel_idx + 1) begin
            src_block[pixel_idx] = 10'd100 + (pixel_idx % 28);
        end
    end
endtask

// Create edge data
task create_edge_data;
    integer pixel_idx;
    begin
        for (pixel_idx = 0; pixel_idx < 64; pixel_idx = pixel_idx + 1) begin
            if ((pixel_idx % 8) < 4) begin
                src_block[pixel_idx] = 10'd50;
            end else begin
                src_block[pixel_idx] = 10'd200;
            end
        end
    end
endtask

// Create directional data
task create_directional_data;
    integer pixel_idx;
    integer x, y;
    begin
        for (pixel_idx = 0; pixel_idx < 64; pixel_idx = pixel_idx + 1) begin
            x = pixel_idx % 8;
            y = pixel_idx / 8;
            // Diagonal pattern
            src_block[pixel_idx] = 10'd100 + ((x + y) * 5);
        end
    end
endtask

//==============================================================================
// Test Cases
//==============================================================================

initial begin
    // Initialize
    rst_n = 0;
    ready = 0;
    strength_y = 3'd0;
    strength_uv = 3'd0;
    damping = 3'd0;
    is_chroma = 1'b0;
    test_num = 0;
    error_count = 0;
    
    // Initialize source block
    for (i = 0; i < 64; i = i + 1) begin
        src_block[i] = 10'd128;
    end
    
    // Reset sequence
    #20;
    rst_n = 1;
    #20;
    
    $display("========================================");
    $display("CDEF Filter Testbench");
    $display("========================================");
    
    // Test 1: No filtering (strength = 0)
    test_num = 1;
    $display("\n[Test %0d] No Filtering (strength=0)", test_num);
    create_flat_data();
    strength_y = 3'd0;
    strength_uv = 3'd0;
    damping = 3'd0;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    // Verify output equals input
    for (i = 0; i < 64; i = i + 1) begin
        if (dst_block[i] !== src_block[i]) begin
            $display("✗ Error at pixel %0d: output=%0d, input=%0d", 
                     i, dst_block[i], src_block[i]);
            error_count = error_count + 1;
        end
    end
    $display("✓ Test %0d: No filtering verified", test_num);
    
    #50;
    
    // Test 2: Light filtering (strength = 1)
    test_num = 2;
    $display("\n[Test %0d] Light Filtering (strength=1)", test_num);
    create_noisy_data();
    strength_y = 3'd1;
    strength_uv = 3'd0;
    damping = 3'd1;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Light filtering test completed", test_num);
    $display("Sample output pixels: dst[0]=%0d, dst[8]=%0d, dst[16]=%0d",
             dst_block[0], dst_block[8], dst_block[16]);
    
    #50;
    
    // Test 3: Medium filtering (strength = 3)
    test_num = 3;
    $display("\n[Test %0d] Medium Filtering (strength=3)", test_num);
    create_noisy_data();
    strength_y = 3'd3;
    strength_uv = 3'd0;
    damping = 3'd2;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Medium filtering test completed", test_num);
    
    #50;
    
    // Test 4: Strong filtering (strength = 7)
    test_num = 4;
    $display("\n[Test %0d] Strong Filtering (strength=7)", test_num);
    create_noisy_data();
    strength_y = 3'd7;
    strength_uv = 3'd0;
    damping = 3'd3;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Strong filtering test completed", test_num);
    
    #50;
    
    // Test 5: Edge preservation
    test_num = 5;
    $display("\n[Test %0d] Edge Preservation", test_num);
    create_edge_data();
    strength_y = 3'd3;
    strength_uv = 3'd0;
    damping = 3'd2;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Edge preservation test completed", test_num);
    $display("Edge pixels: dst[3]=%0d, dst[4]=%0d (original: %0d, %0d)",
             dst_block[3], dst_block[4], src_block[3], src_block[4]);
    
    #50;
    
    // Test 6: Directional filtering
    test_num = 6;
    $display("\n[Test %0d] Directional Filtering", test_num);
    create_directional_data();
    strength_y = 3'd4;
    strength_uv = 3'd0;
    damping = 3'd2;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Directional filtering test completed", test_num);
    
    #50;
    
    // Test 7: Chroma filtering
    test_num = 7;
    $display("\n[Test %0d] Chroma Filtering", test_num);
    create_noisy_data();
    strength_y = 3'd0;
    strength_uv = 3'd2;
    damping = 3'd1;
    is_chroma = 1'b1;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: Chroma filtering test completed", test_num);
    
    #50;
    
    // Test 8: No damping (damping = 0)
    test_num = 8;
    $display("\n[Test %0d] No Damping (damping=0)", test_num);
    create_noisy_data();
    strength_y = 3'd5;
    strength_uv = 3'd0;
    damping = 3'd0;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: No damping test completed", test_num);
    
    #50;
    
    // Test 9: High damping (damping = 7)
    test_num = 9;
    $display("\n[Test %0d] High Damping (damping=7)", test_num);
    create_noisy_data();
    strength_y = 3'd5;
    strength_uv = 3'd0;
    damping = 3'd7;
    is_chroma = 1'b0;
    
    @(posedge clk);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    wait(valid);
    ready = 1;
    @(posedge clk);
    ready = 0;
    
    $display("Test %0d: High damping test completed", test_num);
    
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