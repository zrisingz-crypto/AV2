module test_simple_itx;
parameter CLK_PERIOD = 10;
reg clk;
reg rst_n;
reg start;
reg signed [15:0] coeffs [0:4095];
wire signed [15:0] pixels [0:4095];
wire valid;
wire done;

initial begin
    clk = 0;
    forever #(CLK_PERIOD/2) clk = ~clk;
end

initial begin
    rst_n = 0;
    #(CLK_PERIOD * 5);
    rst_n = 1;
end

av2_inverse_transform_simple u_itx (
    .clk(clk),
    .rst_n(rst_n),
    .coeffs(coeffs),
    .num_coeffs(16'd256),
    .tx_width(6'd16),
    .tx_height(6'd16),
    .tx_type(4'd0),
    .start(start),
    .pixels(pixels),
    .valid(valid),
    .ready(1'b1),
    .done(done)
);

initial begin
    integer i;
    
    // Initialize coeffs
    for (i = 0; i < 4096; i = i + 1)
        coeffs[i] = 16'sd0;
    
    start = 0;
    
    @(posedge rst_n);
    @(posedge clk);
    
    // Set some coefficient values
    coeffs[0] = 16'sd1;
    coeffs[1] = 16'sd2;
    coeffs[2] = 16'sd3;
    
    $display("Before start: coeffs[0]=%0d, coeffs[1]=%0d", coeffs[0], coeffs[1]);
    
    @(posedge clk);
    start = 1;
    @(posedge clk);
    start = 0;
    
    wait(done);
    @(posedge clk);
    
    $display("After done: pixels[0]=%0d, pixels[1]=%0d", pixels[0], pixels[1]);
    
    $finish;
end

endmodule
