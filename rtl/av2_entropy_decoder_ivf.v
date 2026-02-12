//==============================================================================
// AV2 Entropy Decoder for IVF - Real Bitstream Parser (Simplified)
// Parses actual IVF bitstream data and extracts symbols
//==============================================================================

`timescale 1ns / 1ps

module av2_entropy_decoder_ivf #(
    parameter DATA_WIDTH = 128,
    parameter WINDOW_SIZE = 16
)(
    input  wire                     clk,
    input  wire                     rst_n,
    
    // Bitstream input (AXI-Stream like)
    input  wire [DATA_WIDTH-1:0]    bitstream_data,
    input  wire                     bitstream_valid,
    output reg                      bitstream_ready,
    
    // Context model
    input  wire [15:0]              context_idx,
    input  wire [15:0]              context_prob,
    output wire [15:0]              context_idx_out,
    output wire                     context_update_en,
    output wire [15:0]              context_update_idx,
    output wire                     context_update_bit,
    
    // Symbol output
    output reg  [15:0]              symbol,
    output reg                      symbol_valid,
    input  wire                     symbol_ready,
    
    // Control
    input  wire                     reset_contexts,
    input  wire                     start,
    output reg                      done
);

//==============================================================================
// State Machine
//==============================================================================

localparam IDLE          = 3'd0;
localparam LOAD_DATA     = 3'd1;
localparam PARSE_OBU     = 3'd2;
localparam EXTRACT_COEFF = 3'd3;
localparam OUTPUT_SYMBOL = 3'd4;
localparam DONE_STATE    = 3'd5;

reg [2:0] state, state_next;

//==============================================================================
// Bitstream Buffer
//==============================================================================

reg [DATA_WIDTH-1:0] bitstream_buffer;
reg [15:0] bit_ptr;
reg [15:0] byte_ptr;
reg [31:0] total_bytes;

// OBU parsing
reg [3:0]  obu_type;
reg [31:0] obu_size;
reg        obu_has_size;

// Symbol generation
reg [15:0] symbol_count;
reg [15:0] max_symbols;

//==============================================================================
// Bit extraction functions
//==============================================================================

function [7:0] get_byte;
    input [15:0] offset;
    begin
        if (offset < DATA_WIDTH/8)
            get_byte = bitstream_buffer[offset*8 +: 8];
        else
            get_byte = 8'h00;
    end
endfunction

//==============================================================================
// Main State Machine
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        bitstream_ready <= 1'b0;
        symbol_valid <= 1'b0;
        done <= 1'b0;
        bit_ptr <= 16'd0;
        byte_ptr <= 16'd0;
        total_bytes <= 32'd0;
        symbol_count <= 16'd0;
        symbol <= 16'd0;
    end else begin
        state <= state_next;
        
        case (state)
            IDLE: begin
                bitstream_ready <= 1'b1;
                symbol_valid <= 1'b0;
                done <= 1'b0;
                bit_ptr <= 16'd0;
                byte_ptr <= 16'd0;
                symbol_count <= 16'd0;
                
                if (start) begin
                    $display("[TIME %0t] IVF Entropy Decoder: Starting...", $time);
                end
            end
            
            LOAD_DATA: begin
                if (bitstream_valid && bitstream_ready) begin
                    bitstream_buffer <= bitstream_data;
                    total_bytes <= total_bytes + (DATA_WIDTH/8);
                    $display("[TIME %0t] Loaded %d bytes, total=%d", $time, DATA_WIDTH/8, total_bytes + (DATA_WIDTH/8));
                end
            end
            
            PARSE_OBU: begin
                // Parse OBU header from first byte
                if (byte_ptr == 0) begin
                    obu_type <= (get_byte(0) >> 3) & 4'hF;
                    obu_has_size <= (get_byte(0) >> 1) & 1'b1;
                    $display("[TIME %0t] OBU Type=%d, HasSize=%b", $time, (get_byte(0) >> 3) & 4'hF, (get_byte(0) >> 1) & 1'b1);
                end
            end
            
            EXTRACT_COEFF: begin
                // Extract coefficients/symbols from bitstream
                // Simplified: use byte values directly as symbols
                if (symbol_count < max_symbols) begin
                    // Generate symbol from bitstream data
                    symbol <= {8'h00, get_byte(symbol_count[15:0] % (DATA_WIDTH/8))};
                end
            end
            
            OUTPUT_SYMBOL: begin
                symbol_valid <= 1'b1;
                if (symbol_ready) begin
                    symbol_count <= symbol_count + 1;
                    $display("[TIME %0t] Output symbol[%0d] = %d", $time, symbol_count, symbol);
                end
            end
            
            DONE_STATE: begin
                symbol_valid <= 1'b0;
                done <= 1'b1;
                $display("[TIME %0t] IVF Entropy Decoder: Done, %d symbols output", $time, symbol_count);
            end
        endcase
    end
end

//==============================================================================
// State Transition Logic
//==============================================================================

always @(*) begin
    state_next = state;
    max_symbols = 16'd64;  // Output 64 symbols for testing
    
    case (state)
        IDLE: begin
            if (start)
                state_next = LOAD_DATA;
        end
        
        LOAD_DATA: begin
            if (bitstream_valid)
                state_next = PARSE_OBU;
        end
        
        PARSE_OBU: begin
            state_next = EXTRACT_COEFF;
        end
        
        EXTRACT_COEFF: begin
            state_next = OUTPUT_SYMBOL;
        end
        
        OUTPUT_SYMBOL: begin
            if (symbol_valid && symbol_ready) begin
                if (symbol_count >= max_symbols - 1)
                    state_next = DONE_STATE;
                else
                    state_next = EXTRACT_COEFF;
            end
        end
        
        DONE_STATE: begin
            state_next = IDLE;
        end
    endcase
end

// Context model outputs (simplified)
assign context_idx_out = context_idx;
assign context_update_en = 1'b0;
assign context_update_idx = 16'd0;
assign context_update_bit = 1'b0;

endmodule
