//==============================================================================
// Real Entropy Decoder Module for AV2
// Based on range coder algorithm from software implementation
//==============================================================================

`timescale 1ns / 1ps

module av2_entropy_decoder_real #(
    parameter DATA_WIDTH = 128,
    parameter WINDOW_SIZE = 32,
    parameter CDF_PROB_TOP = 32768
)(
    input  wire                      clk,
    input  wire                      rst_n,
    
    // Bitstream input
    input  wire [DATA_WIDTH-1:0]    bitstream_data,
    input  wire                      bitstream_valid,
    output reg                       bitstream_ready,
    
    // Context model interface
    input  wire [15:0]              context_idx,
    input  wire [15:0]              context_prob,
    output reg  [15:0]              context_idx_out,
    output reg                       context_update_en,
    output reg  [15:0]              context_update_idx,
    output reg                       context_update_bit,
    input  wire                      reset_contexts,
    
    // Symbol output
    output reg  [15:0]              symbol,
    output reg                       symbol_valid,
    input  wire                      symbol_ready,
    
    // Control
    input  wire                      start,
    output wire                      done
);

// Range decoder state (based on od_ec_dec in entdec.c)
reg [WINDOW_SIZE-1:0] dif;           // Decoded difference window
reg [15:0]             rng;            // Range
reg [15:0]             cnt;            // Bit count

// Bitstream buffer
reg [DATA_WIDTH*8-1:0]  bitstream_buffer;  // Internal buffer
reg [DATA_WIDTH-1:0]    bitstream_reg;
reg [15:0]             bit_ptr;        // Pointer in bitstream
reg [15:0]             bit_count;      // Total bits processed

// State machine
localparam IDLE       = 2'd0;
localparam DECODING   = 2'd1;
localparam REFILL     = 2'd2;
localparam DONE       = 2'd3;

reg [1:0] state, state_next;

// Decoding registers
reg [4:0]  symbol_idx;       // Index for decoding multiple symbols
reg [4:0]  total_symbols;    // Total symbols to decode

// Temporary calculation registers
reg [WINDOW_SIZE-1:0] vw;     // Value window
reg [15:0]             r;       // Range
reg [15:0]             c;       // Cumulative distribution
reg [15:0]             u, v;    // Intermediate values
reg [15:0]             ret;     // Decoded symbol

// Initialization constants
initial begin
    dif = (1'b1 << (WINDOW_SIZE - 1)) - 1;  // Initialize dif
    rng = 16'h8000;                             // Initialize rng
    cnt = 16'hFFF1;                             // Initialize cnt
end

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        bitstream_ready <= 1'b0;
        symbol <= 16'd0;
        symbol_valid <= 1'b0;
        symbol_idx <= 5'd0;
        total_symbols <= 5'd0;
        bit_ptr <= 16'd0;
        bit_count <= 16'd0;
        bitstream_buffer <= {(DATA_WIDTH*8){1'b0}};
        bitstream_reg <= {(DATA_WIDTH){1'b0}};
        context_idx_out <= 16'd0;
        context_update_en <= 1'b0;
        context_update_idx <= 16'd0;
        context_update_bit <= 1'b0;
    end else begin
        state <= state_next;
        
        case (state)
            IDLE: begin
                bitstream_ready <= 1'b1;
                symbol_valid <= 1'b0;
                
                if (start) begin
                    // Initialize decoder state
                    dif <= (1'b1 << (WINDOW_SIZE - 1)) - 1;
                    rng <= 16'h8000;
                    cnt <= -15;  // Initial count for refill
                    bit_ptr <= 16'd0;
                    bit_count <= 16'd0;
                    symbol_idx <= 5'd0;
                    
                    // Assume 8 symbols to decode for now
                    total_symbols <= 5'd8;
                    
                    state <= REFILL;
                end
            end
            
            REFILL: begin
                // Refill the difference window with bits from bitstream
                if (bitstream_valid && bitstream_ready) begin
                    // Load new bitstream data
                    bitstream_reg <= bitstream_data;
                    bitstream_buffer <= {bitstream_buffer[DATA_WIDTH*8-1-DATA_WIDTH:0], 
                                         bitstream_data};
                    bit_count <= bit_count + DATA_WIDTH;
                    
                    // Refill dif (simplified version)
                    if (cnt >= 16'd0) begin
                        dif <= dif << 8;
                        cnt <= cnt - 8'd8;
                    end
                    
                    // Always transition to DECODING after refilling
                    state <= DECODING;
                end
            end
            
            DECODING: begin
                symbol_valid <= 1'b0;
                
                if (symbol_idx < total_symbols) begin
                    // Decode symbol using context model
                    // Simplified CDF-based decoding
                    
                    // Calculate normalized values
                    r = rng;
                    c = dif >> (WINDOW_SIZE - 16);  // Compare value
                    
                    // Simple binary decision based on context_prob
                    if (c >= context_prob) begin
                        // Symbol 1
                        ret = 16'd1;
                        u = context_prob;
                        v = r - context_prob;
                        dif = dif - (context_prob << (WINDOW_SIZE - 16));
                    end else begin
                        // Symbol 0
                        ret = 16'd0;
                        u = 16'd0;
                        v = context_prob;
                    end
                    
                    // Update range
                    rng = v;
                    
                    // Output symbol
                    if (symbol_ready || !symbol_valid) begin
                        symbol <= ret;
                        symbol_valid <= 1'b1;
                        symbol_idx <= symbol_idx + 1;
                        
                        // Update context
                        context_idx_out <= context_idx;
                        context_update_en <= 1'b1;
                        context_update_idx <= context_idx;
                        context_update_bit <= ret[0];
                    end
                    
                    // Normalize (refill)
                    if (rng < 16'h8000) begin
                        // Need to refill range
                        dif <= dif << 1;
                        rng <= rng << 1;
                        cnt <= cnt + 1;
                        
                        if (cnt >= 16'd0) begin
                            state <= REFILL;
                        end
                    end
                end else begin
                    // Done decoding
                    state <= DONE;
                end
            end
            
            DONE: begin
                symbol_valid <= 1'b0;
                if (symbol_ready) begin
                    state <= IDLE;
                end
            end
            
            default: begin
                state <= IDLE;
            end
        endcase
        
        // Handle context reset
        if (reset_contexts) begin
            context_update_en <= 1'b0;
        end
    end
end

// State transition logic
always @(*) begin
    state_next = state;
    
    case (state)
        IDLE: begin
            if (start)
                state_next = REFILL;
        end
        
        REFILL: begin
            if (bitstream_valid && bitstream_ready && cnt >= 16'd0)
                state_next = DECODING;
        end
        
        DECODING: begin
            if (symbol_idx >= total_symbols)
                state_next = DONE;
            else if (rng < 16'h8000 && cnt >= 16'd0)
                state_next = REFILL;
        end
        
        DONE: begin
            if (symbol_ready)
                state_next = IDLE;
        end
        
        default: begin
            state_next = IDLE;
        end
    endcase
end

assign done = (state == DONE);

endmodule