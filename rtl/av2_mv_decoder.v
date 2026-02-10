//==============================================================================
// Motion Vector Decoder Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_mv_decoder (
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [15:0]             context_idx,
    input  wire [15:0]             context_prob,
    input  wire [15:0]             decoded_symbol,
    input  wire                    symbol_valid,
    output reg                     symbol_ready,
    output reg  signed [15:0]      mv_x,
    output reg  signed [15:0]      mv_y,
    output reg                     mv_valid,
    input  wire                    mv_ready,
    input  wire                    start,
    output wire                    done
);

// Simple stub: return zero MV
reg busy;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        symbol_ready <= 1'b0;
        mv_x <= 16'sd0;
        mv_y <= 16'sd0;
        mv_valid <= 1'b0;
        busy <= 1'b0;
    end else begin
        if (start && !busy) begin
            busy <= 1'b1;
            symbol_ready <= 1'b1;
        end else if (busy) begin
            if (symbol_valid) begin
                symbol_ready <= 1'b0;
                mv_x <= 16'sd0;
                mv_y <= 16'sd0;
                mv_valid <= 1'b1;
            end else if (mv_ready) begin
                mv_valid <= 1'b0;
                busy <= 1'b0;
            end
        end
    end
end

assign done = !busy;

endmodule