//==============================================================================
// Frame Header Parser Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_frame_header_parser (
    input  wire                    clk,
    input  wire                    rst_n,
    
    input  wire                    obu_valid,
    input  wire [127:0]            obu_data,
    
    output wire [1:0]              frame_type,
    output wire [15:0]             frame_width,
    output wire [15:0]             frame_height,
    output wire [7:0]              qindex,
    output wire                    header_valid,
    input  wire                    header_ready
);

// Simple stub: return default frame parameters
reg [1:0] frame_type_r;
reg [15:0] frame_width_r;
reg [15:0] frame_height_r;
reg [7:0] qindex_r;
reg header_valid_r;
reg parsed_flag;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        frame_type_r <= 2'd0;
        frame_width_r <= 16'd64;
        frame_height_r <= 16'd64;
        qindex_r <= 8'd32;
        header_valid_r <= 1'b0;
        parsed_flag <= 1'b0;
    end else begin
        // Parse header only once on valid OBU
        if (obu_valid && !parsed_flag && !header_valid_r) begin
            frame_type_r <= 2'd0;  // KEY_FRAME
            frame_width_r <= 16'd64;
            frame_height_r <= 16'd64;
            qindex_r <= 8'd32;
            header_valid_r <= 1'b1;
            parsed_flag <= 1'b1;
        end else if (header_ready && header_valid_r) begin
            header_valid_r <= 1'b0;
        end
        
        // Reset parse flag when OBU is no longer valid
        if (!obu_valid) begin
            parsed_flag <= 1'b0;
        end
    end
end

assign frame_type = frame_type_r;
assign frame_width = frame_width_r;
assign frame_height = frame_height_r;
assign qindex = qindex_r;
assign header_valid = header_valid_r;

endmodule