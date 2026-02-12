//==============================================================================
// IVF Bitstream Data Module
// Source: crowd_64x64p30f32_av2.ivf
// Resolution: 64x64
// Frames extracted: 2
//==============================================================================

`ifndef IVF_BITSTREAM_DATA_V
`define IVF_BITSTREAM_DATA_V

module ivf_bitstream_data #(
    parameter FRAME_0_SIZE = 3392,
    parameter FRAME_1_SIZE = 47,
    parameter IVF_NUM_FRAMES = 2
)(
    output reg [7:0] frame_0_data [0:FRAME_0_SIZE-1],
    output reg [7:0] frame_1_data [0:FRAME_1_SIZE-1]
);

initial begin
    $display("[INFO] Loading IVF bitstream data...");
    $readmemh("ivf_bitstream_data_frame_0.hex", frame_0_data);
    $readmemh("ivf_bitstream_data_frame_1.hex", frame_1_data);
    $display("[INFO] IVF bitstream data loaded successfully");
    $display("[INFO] Frame 0: %0d bytes", FRAME_0_SIZE);
    $display("[INFO] Frame 1: %0d bytes", FRAME_1_SIZE);
end

endmodule

`endif
