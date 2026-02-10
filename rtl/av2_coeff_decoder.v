//==============================================================================
// Coefficient Decoder Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_coeff_decoder #(
    parameter MAX_COEFFS = 4096
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [15:0]             context_idx,
    input  wire [15:0]             context_prob,
    input  wire [15:0]             decoded_symbol,
    input  wire                    symbol_valid,
    output reg                     symbol_ready,
    output reg  signed [15:0]      coeffs[0:MAX_COEFFS-1],
    output reg  [15:0]             num_coeffs,
    output reg                     coeffs_valid,
    input  wire                    coeffs_ready,
    input  wire [5:0]              tx_size,
    input  wire                    start,
    output wire                    done
);

// Simple stub: Simulate coefficient decoding with fixed latency
reg [7:0] coeff_count;
reg [7:0] cycle_count;
reg busy;

integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        symbol_ready <= 1'b0;
        coeffs_valid <= 1'b0;
        num_coeffs <= 16'd0;
        coeff_count <= 8'd0;
        cycle_count <= 8'd0;
        busy <= 1'b0;
        for (i = 0; i < MAX_COEFFS; i = i + 1) begin
            coeffs[i] <= 16'sd0;
        end
    end else begin
        if (start && !busy) begin
            busy <= 1'b1;
            coeff_count <= 8'd0;
            cycle_count <= 8'd0;
            symbol_ready <= 1'b1;
        end else if (busy) begin
            // Fixed latency approach - complete after 5 cycles regardless of input
            cycle_count <= cycle_count + 1;
            
            if (cycle_count == 8'd5) begin
                coeffs_valid <= 1'b1;
                symbol_ready <= 1'b0;
            end else if (cycle_count > 8'd5) begin
                num_coeffs <= 16'd16;
                if (coeffs_ready) begin
                    coeffs_valid <= 1'b0;
                    busy <= 1'b0;
                    coeff_count <= 8'd0;
                    cycle_count <= 8'd0;
                end
            end
        end
    end
end

assign done = !busy;

endmodule