// Simple debug module to trace state machine
module debug_decoder;
initial begin
    $dumpfile("decoder.vcd");
    $dumpvars(0, tb_tile_decoder_real_fixed);
end
endmodule
