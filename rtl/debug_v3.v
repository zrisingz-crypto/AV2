// Debug module to trace x values
module debug_trace;
initial begin
    $dumpfile("debug_v3.vcd");
    $dumpvars(0, tb_tile_decoder_v3);
end
endmodule
