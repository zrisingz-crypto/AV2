//==============================================================================
// Testbench for Tile Decoder V4
//==============================================================================

`timescale 1ns / 1ps

module tb_tile_decoder_v4;

parameter CLK_PERIOD = 10;
parameter MAX_WIDTH = 64;
parameter MAX_HEIGHT = 64;
parameter TIMEOUT_CYCLES = 10000;

reg clk;
reg rst_n;
reg start;
reg [15:0] frame_width;
reg [15:0] frame_height;
reg [7:0] qindex;
reg [1:0] frame_type;
reg [127:0] tile_data;
reg tile_valid;
wire tile_ready;
wire [31:0] ref_read_addr;
reg [9:0] ref_pixel_data;
wire ref_read_en;
wire [127:0] recon_data;
wire [31:0] recon_addr;
wire recon_wr_en;
wire tile_done;

integer cycle_count;
integer output_file;

// DUT
av2_tile_decoder_v4 #(
    .MAX_WIDTH(MAX_WIDTH),
    .MAX_HEIGHT(MAX_HEIGHT)
) dut (
    .clk(clk),
    .rst_n(rst_n),
    .start(start),
    .frame_width(frame_width),
    .frame_height(frame_height),
    .qindex(qindex),
    .frame_type(frame_type),
    .tile_data(tile_data),
    .tile_valid(tile_valid),
    .tile_ready(tile_ready),
    .ref_read_addr(ref_read_addr),
    .ref_pixel_data(ref_pixel_data),
    .ref_read_en(ref_read_en),
    .recon_data(recon_data),
    .recon_addr(recon_addr),
    .recon_wr_en(recon_wr_en),
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

// Cycle counter
always @(posedge clk) begin
    if (!rst_n)
        cycle_count <= 0;
    else
        cycle_count <= cycle_count + 1;
end

// Test
initial begin
    output_file = $fopen("rtl_v4_output.txt", "w");
    
    start = 0;
    frame_width = 64;
    frame_height = 64;
    qindex = 128;
    frame_type = 0;
    tile_data = 0;
    tile_valid = 0;
    ref_pixel_data = 128;
    
    @(posedge rst_n);
    @(posedge clk);
    @(posedge clk);
    
    $display("\n========================================");
    $display("Tile Decoder V4 Test");
    $display("========================================\n");
    
    @(posedge clk);
    tile_valid <= 1;
    tile_data <= {32'h12345678, 32'h9ABCDEF0, 32'h13579BDF, 32'h2468ACE0};
    
    @(posedge clk);
    start <= 1;
    @(posedge clk);
    start <= 0;
    
    repeat(20) begin
        @(posedge clk);
        if (tile_ready) begin
            tile_valid <= 1;
            tile_data <= tile_data + 1;
        end
    end
    
    tile_valid <= 0;
    
    wait(tile_done);
    @(posedge clk);
    $display("[TIME %0t] Decode complete! Total cycles: %0d", $time, cycle_count);
    
    // Print some recon_frame values
    $display("\nRecon frame values [0:63]:");
    for (integer i = 0; i < 64; i = i + 1) begin
        $write("%0d ", dut.recon_frame[i]);
        if (i % 16 == 15) $display("");
    end
    
    @(posedge clk);
    @(posedge clk);
    
    $fclose(output_file);
    
    $display("\n========================================");
    $display("Test Done!");
    $display("========================================\n");
    
    $finish;
end

// Monitor output writes
always @(posedge clk) begin
    if (recon_wr_en) begin
        $display("[TIME %0t] Write: addr=%0d, data[0]=%0d", $time, recon_addr, recon_data[7:0]);
        
        for (integer j = 0; j < 16; j = j + 1) begin
            reg [7:0] pixel;
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
            $fwrite(output_file, "pixel[%0d][%0d] = %0d\n", recon_addr, j, pixel);
        end
    end
end

// Timeout check
initial begin
    #(CLK_PERIOD * TIMEOUT_CYCLES);
    $display("[ERROR] Timeout after %0d cycles!", TIMEOUT_CYCLES);
    $fclose(output_file);
    $finish;
end

endmodule
