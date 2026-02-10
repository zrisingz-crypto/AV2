//==============================================================================
// Entropy Decoder Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_entropy_decoder #(
    parameter DATA_WIDTH = 128
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [DATA_WIDTH-1:0]   bitstream_data,
    input  wire                    bitstream_valid,
    output reg                     bitstream_ready,
    input  wire [15:0]             context_idx,
    input  wire [15:0]             context_prob,
    output reg  [15:0]             symbol,
    output reg                     symbol_valid,
    input  wire                    symbol_ready,
    input  wire                    start,
    output wire                    done
);

// Simple stub: return fixed symbols
reg [3:0] decode_count;
reg busy;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        bitstream_ready <= 1'b1;
        symbol <= 16'd0;
        symbol_valid <= 1'b0;
        decode_count <= 4'd0;
        busy <= 1'b0;
    end else begin
        if (start && !busy) begin
            busy <= 1'b1;
            decode_count <= 4'd0;
        end else if (busy) begin
            if (decode_count < 4'd8 && symbol_ready) begin
                symbol <= {4'd0, decode_count, 8'd0}; // Simple symbol pattern
                symbol_valid <= 1'b1;
                decode_count <= decode_count + 1;
            end else if (decode_count >= 4'd8) begin
                symbol_valid <= 1'b0;
                busy <= 1'b0;
            end else if (!symbol_valid) begin
                symbol_valid <= 1'b1;
            end else if (symbol_ready) begin
                symbol_valid <= 1'b0;
            end
        end
    end
end

assign done = !busy && (decode_count >= 4'd8);

endmodule