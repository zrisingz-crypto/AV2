//==============================================================================
// Context Model Module (Stub for Simulation)
//==============================================================================

`timescale 1ns / 1ps

module av2_context_model #(
    parameter NUM_CONTEXTS = 1024
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [15:0]             context_idx,
    output reg  [15:0]             context_prob,
    input  wire                    update_en,
    input  wire [15:0]             update_idx,
    input  wire                    update_bit,
    input  wire                    reset_contexts
);

// Simple context model with fixed probabilities
reg [15:0] context_probs[0:NUM_CONTEXTS-1];

integer i;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (i = 0; i < NUM_CONTEXTS; i = i + 1) begin
            context_probs[i] <= 16'd128; // Default probability (50%)
        end
        context_prob <= 16'd128;
    end else begin
        // Update probability on valid update
        if (update_en && update_idx < NUM_CONTEXTS) begin
            if (update_bit)
                context_probs[update_idx] <= context_probs[update_idx] + 16'd1;
            else
                context_probs[update_idx] <= context_probs[update_idx] - 16'd1;
        end
        
        // Reset all contexts
        if (reset_contexts) begin
            for (i = 0; i < NUM_CONTEXTS; i = i + 1) begin
                context_probs[i] <= 16'd128;
            end
        end
        
        // Output current context probability
        if (context_idx < NUM_CONTEXTS)
            context_prob <= context_probs[context_idx];
        else
            context_prob <= 16'd128;
    end
end

endmodule