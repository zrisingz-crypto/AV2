//==============================================================================
// OBU Parser Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_obu_parser #(
    parameter DATA_WIDTH = 128
)(
    input  wire                    clk,
    input  wire                    rst_n,
    
    // AXI4-Stream Input
    input  wire [DATA_WIDTH-1:0]   s_axis_tdata,
    input  wire                    s_axis_tvalid,
    output wire                    s_axis_tready,
    input  wire                    s_axis_tlast,
    
    // OBU Output
    output wire [3:0]              obu_type,
    output wire [31:0]             obu_size,
    output wire                    obu_valid,
    input  wire                    obu_ready
);

// Simple stub: detect OBU headers in incoming stream
reg [3:0] obu_type_r;
reg [31:0] obu_size_r;
reg obu_valid_r;
reg s_axis_tready_r;
reg obu_parsed;  // Flag to indicate OBU has been parsed

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        obu_type_r <= 4'd0;
        obu_size_r <= 32'd0;
        obu_valid_r <= 1'b0;
        s_axis_tready_r <= 1'b1;
        obu_parsed <= 1'b0;
    end else begin
        // Accept input when ready and not already parsed
        if (s_axis_tvalid && s_axis_tready_r && !obu_parsed) begin
            // Extract OBU type from first byte (bits 7-4 are obu_type)
            obu_type_r <= s_axis_tdata[7:4];
            obu_size_r <= s_axis_tdata[127:96]; // Simple size extraction
            obu_valid_r <= 1'b1;
            obu_parsed <= 1'b1;
            s_axis_tready_r <= 1'b0; // Not ready until consumed
        end else if (obu_valid_r && obu_ready) begin
            // OBU consumed, clear valid and reset for next
            obu_valid_r <= 1'b0;
            obu_parsed <= 1'b0;
            s_axis_tready_r <= 1'b1; // Ready for next input
        end
    end
end

assign obu_type = obu_type_r;
assign obu_size = obu_size_r;
assign obu_valid = obu_valid_r;
assign s_axis_tready = s_axis_tready_r;

endmodule