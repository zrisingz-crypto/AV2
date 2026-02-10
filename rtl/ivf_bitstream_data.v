//==============================================================================
// IVF Bitstream Data
// Source: crowd_64x64p30f32_av2.ivf
// Resolution: 64x64
// Frames extracted: 2
//==============================================================================

`define IVF_WIDTH 64
`define IVF_HEIGHT 64
`define IVF_NUM_FRAMES 2
`define IVF_HEADER_SIZE 12

// Frame 0: 3392 bytes
localparam FRAME_0_SIZE = 3392;
reg [7:0] frame_0_data [0:3392-1];
initial $readmemh("ivf_bitstream_data_frame_0.hex", frame_0_data);

// Frame 1: 47 bytes
localparam FRAME_1_SIZE = 47;
reg [7:0] frame_1_data [0:47-1];
initial $readmemh("ivf_bitstream_data_frame_1.hex", frame_1_data);


// Frame sizes array
reg [31:0] frame_sizes [0:1];
initial begin
    frame_sizes[0] = FRAME_0_SIZE;
    frame_sizes[1] = FRAME_1_SIZE;
end

// Total bitstream size
localparam BITSTREAM_SIZE = 3439;
