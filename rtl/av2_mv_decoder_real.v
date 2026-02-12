//==============================================================================
// Motion Vector Decoder Module (Real Implementation)
// Decodes motion vectors from entropy decoder output
//==============================================================================

`timescale 1ns / 1ps

module av2_mv_decoder_real (
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

//==============================================================================
// Parameters
//==============================================================================

localparam STATE_IDLE       = 3'd0;
localparam STATE_DECODE_X   = 3'd1;
localparam STATE_DECODE_Y   = 3'd2;
localparam STATE_SIGN_X     = 3'd3;
localparam STATE_SIGN_Y     = 3'd4;
localparam STATE_OUTPUT     = 3'd5;
localparam STATE_DONE       = 3'd6;

// Motion vector component decoding parameters
localparam MV_MAX_BITS    = 6'd12;  // Maximum bits for MV magnitude
localparam MV_SIGN_BIT    = 1;      // Sign bit position

//==============================================================================
// State Register
//==============================================================================

reg [2:0] state;
reg [2:0] state_next;

//==============================================================================
// Working Registers
//==============================================================================

reg [MV_MAX_BITS-1:0] mv_x_mag;
reg [MV_MAX_BITS-1:0] mv_y_mag;
reg                    mv_x_sign;
reg                    mv_y_sign;
reg [3:0]             bit_count;
reg                    busy;

//==============================================================================
// State Machine - Sequential Logic
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= STATE_IDLE;
        symbol_ready <= 1'b0;
        mv_x <= 16'sd0;
        mv_y <= 16'sd0;
        mv_valid <= 1'b0;
        busy <= 1'b0;
        mv_x_mag <= {MV_MAX_BITS{1'b0}};
        mv_y_mag <= {MV_MAX_BITS{1'b0}};
        mv_x_sign <= 1'b0;
        mv_y_sign <= 1'b0;
        bit_count <= 4'd0;
    end else begin
        state <= state_next;
        
        case (state)
            STATE_IDLE: begin
                if (start) begin
                    busy <= 1'b1;
                    symbol_ready <= 1'b1;
                    mv_x_mag <= {MV_MAX_BITS{1'b0}};
                    mv_y_mag <= {MV_MAX_BITS{1'b0}};
                    mv_x_sign <= 1'b0;
                    mv_y_sign <= 1'b0;
                    bit_count <= 4'd0;
                end
            end
            
            STATE_DECODE_X: begin
                if (symbol_valid) begin
                    // Decode MV magnitude bits using unary-like coding
                    mv_x_mag <= mv_x_mag | (decoded_symbol[bit_count] << bit_count);
                    bit_count <= bit_count + 1;
                    
                    // Check if we've decoded all bits (simplified)
                    if (bit_count >= 4'd11) begin
                        state_next <= STATE_SIGN_X;
                    end else if (decoded_symbol == 16'd0) begin
                        // Zero indicates end of magnitude coding
                        state_next <= STATE_SIGN_X;
                    end
                end
            end
            
            STATE_DECODE_Y: begin
                if (symbol_valid) begin
                    // Decode MV magnitude bits for Y component
                    mv_y_mag <= mv_y_mag | (decoded_symbol[bit_count] << bit_count);
                    bit_count <= bit_count + 1;
                    
                    // Check if we've decoded all bits
                    if (bit_count >= 4'd11) begin
                        state_next <= STATE_SIGN_Y;
                    end else if (decoded_symbol == 16'd0) begin
                        state_next <= STATE_SIGN_Y;
                    end
                end
            end
            
            STATE_SIGN_X: begin
                if (symbol_valid) begin
                    mv_x_sign <= decoded_symbol[0];
                    mv_x <= mv_x_sign ? -({{4'd0}, mv_x_mag}) : {4'd0, mv_x_mag};
                end
            end
            
            STATE_SIGN_Y: begin
                if (symbol_valid) begin
                    mv_y_sign <= decoded_symbol[0];
                    mv_y <= mv_y_sign ? -({{4'd0}, mv_y_mag}) : {4'd0, mv_y_mag};
                end
            end
            
            STATE_OUTPUT: begin
                mv_valid <= 1'b1;
                if (mv_ready) begin
                    mv_valid <= 1'b0;
                    busy <= 1'b0;
                end
            end
            
            STATE_DONE: begin
                busy <= 1'b0;
            end
        endcase
    end
end

//==============================================================================
// State Machine - Combinational Logic
//==============================================================================

always @(*) begin
    state_next = state;
    
    case (state)
        STATE_IDLE: begin
            if (start && !busy) begin
                state_next = STATE_DECODE_X;
            end
        end
        
        STATE_DECODE_X: begin
            if (bit_count >= 4'd11) begin
                state_next = STATE_SIGN_X;
            end else if (symbol_valid && decoded_symbol == 16'd0) begin
                state_next = STATE_SIGN_X;
            end else begin
                state_next = STATE_DECODE_X;
            end
        end
        
        STATE_SIGN_X: begin
            if (symbol_valid) begin
                state_next = STATE_DECODE_Y;
            end
        end
        
        STATE_DECODE_Y: begin
            if (bit_count >= 4'd11) begin
                state_next = STATE_SIGN_Y;
            end else if (symbol_valid && decoded_symbol == 16'd0) begin
                state_next = STATE_SIGN_Y;
            end else begin
                state_next = STATE_DECODE_Y;
            end
        end
        
        STATE_SIGN_Y: begin
            if (symbol_valid) begin
                state_next = STATE_OUTPUT;
            end
        end
        
        STATE_OUTPUT: begin
            if (mv_valid && mv_ready) begin
                state_next = STATE_DONE;
            end else begin
                state_next = STATE_OUTPUT;
            end
        end
        
        STATE_DONE: begin
            state_next = STATE_IDLE;
        end
        
        default: begin
            state_next = STATE_IDLE;
        end
    endcase
end

//==============================================================================
// Symbol Ready Control
//==============================================================================

always @(*) begin
    case (state)
        STATE_DECODE_X: begin
            symbol_ready = (bit_count < 4'd12) && (bit_count < 4'd11 || (bit_count == 4'd0));
        end
        STATE_DECODE_Y: begin
            symbol_ready = (bit_count < 4'd12) && (bit_count < 4'd11 || (bit_count == 4'd0));
        end
        STATE_SIGN_X: begin
            symbol_ready = 1'b1;
        end
        STATE_SIGN_Y: begin
            symbol_ready = 1'b1;
        end
        default: begin
            symbol_ready = 1'b0;
        end
    endcase
end

//==============================================================================
// Done Signal
//==============================================================================

assign done = !busy;

//==============================================================================
// Debug Output
//==============================================================================

always @(posedge clk) begin
    if (state == STATE_OUTPUT && mv_valid) begin
        $display("[MV_DEC] Decoded MV: mv_x=%0d, mv_y=%0d", mv_x, mv_y);
    end
end

endmodule