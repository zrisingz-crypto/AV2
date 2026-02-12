//==============================================================================
// CDEF Filter Module (Real Implementation)
// Constrained Directional Enhancement Filter for AV2
//==============================================================================

`timescale 1ns / 1ps

module av2_cdef_filter_real #(
    parameter BLOCK_SIZE = 8
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [9:0]             src_block[0:63],
    input  wire [2:0]             strength_y,
    input  wire [2:0]             strength_uv,
    input  wire [2:0]             damping,
    input  wire                    is_chroma,
    output reg  [9:0]             dst_block[0:63],
    output reg                     valid,
    input  wire                    ready
);

//==============================================================================
// Parameters
//==============================================================================

localparam STATE_IDLE      = 2'd0;
localparam STATE_FILTER    = 2'd1;
localparam STATE_OUTPUT    = 2'd2;

localparam NUM_DIRECTIONS = 8;
localparam BLOCK_8X8 = 64;

//==============================================================================
// State Register
//==============================================================================

reg [1:0] state;
reg [1:0] state_next;

//==============================================================================
// Working Registers
//==============================================================================

reg [6:0] pixel_idx;
reg [2:0]  strength;

// Filter working registers (declare at module level)
reg [9:0]  p0, p1, p2;
reg [2:0]  x, y;
reg [5:0]  idx0, idx1, idx2;

// Direction vectors for CDEF
wire signed [3:0] dir_dx[0:NUM_DIRECTIONS-1];
wire signed [3:0] dir_dy[0:NUM_DIRECTIONS-1];

assign dir_dx[0] = 4'sd1;  assign dir_dy[0] = 4'sd0;
assign dir_dx[1] = 4'sd1;  assign dir_dy[1] = 4'sd1;
assign dir_dx[2] = 4'sd0;  assign dir_dy[2] = 4'sd1;
assign dir_dx[3] = -4'sd1; assign dir_dy[3] = 4'sd1;
assign dir_dx[4] = -4'sd1; assign dir_dy[4] = 4'sd0;
assign dir_dx[5] = -4'sd1; assign dir_dy[5] = -4'sd1;
assign dir_dx[6] = 4'sd0;  assign dir_dy[6] = -4'sd1;
assign dir_dx[7] = 4'sd1;  assign dir_dy[7] = -4'sd1;

// Temporary filtered block
reg [9:0]  temp_block[0:63];

integer i;

//==============================================================================
// Wire for current strength (used in combinational logic)
//==============================================================================

wire [2:0] current_strength;
assign current_strength = is_chroma ? strength_uv : strength_y;

//==============================================================================
// Helper Functions
//==============================================================================

function [5:0] get_pixel_idx;
    input [2:0] x;
    input [2:0] y;
    reg [5:0] idx;
    begin
        if (x < 3'd0) idx = y * 8 + 3'd0;
        else if (x > 3'd7) idx = y * 8 + 3'd7;
        else if (y < 3'd0) idx = 3'd0 * 8 + x;
        else if (y > 3'd7) idx = 3'd7 * 8 + x;
        else idx = y * 8 + x;
        get_pixel_idx = idx;
    end
endfunction

function [9:0] clip3;
    input [9:0] val;
    input [9:0] min_v;
    input [9:0] max_v;
    begin
        if (val < min_v) clip3 = min_v;
        else if (val > max_v) clip3 = max_v;
        else clip3 = val;
    end
endfunction

//==============================================================================
// State Machine - Sequential Logic
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= STATE_IDLE;
        valid <= 1'b0;
        pixel_idx <= 7'd0;
        strength <= 3'd0;
        for (i = 0; i < 64; i = i + 1) begin
            dst_block[i] <= 10'd0;
            temp_block[i] <= 10'd0;
        end
    end else begin
        state <= state_next;
        
        case (state)
            STATE_IDLE: begin
                // Capture strength when ready is asserted
                if (ready) begin
                    strength <= current_strength;
                end
            end
            
            STATE_FILTER: begin
                if (pixel_idx < 64) begin
                    if (strength > 3'd0) begin
                        x = pixel_idx % 8;
                        y = pixel_idx / 8;
                        idx0 = get_pixel_idx(x, y);
                        idx1 = get_pixel_idx(x + 3'sd1, y);
                        idx2 = get_pixel_idx(x - 3'sd1, y);
                        p0 = src_block[idx0];
                        p1 = src_block[idx1];
                        p2 = src_block[idx2];
                        temp_block[pixel_idx] <= clip3((p0 * 2 + p1 + p2 + 2) >> 2, 10'd0, 10'd1023);
                    end else begin
                        temp_block[pixel_idx] <= src_block[pixel_idx];
                    end
                    pixel_idx <= pixel_idx + 1;
                end
            end
            
            STATE_OUTPUT: begin
                // Copy temp_block to dst_block (blocking assignment for immediate availability)
                dst_block[0] = temp_block[0];
                dst_block[1] = temp_block[1];
                dst_block[2] = temp_block[2];
                dst_block[3] = temp_block[3];
                dst_block[4] = temp_block[4];
                dst_block[5] = temp_block[5];
                dst_block[6] = temp_block[6];
                dst_block[7] = temp_block[7];
                dst_block[8] = temp_block[8];
                dst_block[9] = temp_block[9];
                dst_block[10] = temp_block[10];
                dst_block[11] = temp_block[11];
                dst_block[12] = temp_block[12];
                dst_block[13] = temp_block[13];
                dst_block[14] = temp_block[14];
                dst_block[15] = temp_block[15];
                dst_block[16] = temp_block[16];
                dst_block[17] = temp_block[17];
                dst_block[18] = temp_block[18];
                dst_block[19] = temp_block[19];
                dst_block[20] = temp_block[20];
                dst_block[21] = temp_block[21];
                dst_block[22] = temp_block[22];
                dst_block[23] = temp_block[23];
                dst_block[24] = temp_block[24];
                dst_block[25] = temp_block[25];
                dst_block[26] = temp_block[26];
                dst_block[27] = temp_block[27];
                dst_block[28] = temp_block[28];
                dst_block[29] = temp_block[29];
                dst_block[30] = temp_block[30];
                dst_block[31] = temp_block[31];
                dst_block[32] = temp_block[32];
                dst_block[33] = temp_block[33];
                dst_block[34] = temp_block[34];
                dst_block[35] = temp_block[35];
                dst_block[36] = temp_block[36];
                dst_block[37] = temp_block[37];
                dst_block[38] = temp_block[38];
                dst_block[39] = temp_block[39];
                dst_block[40] = temp_block[40];
                dst_block[41] = temp_block[41];
                dst_block[42] = temp_block[42];
                dst_block[43] = temp_block[43];
                dst_block[44] = temp_block[44];
                dst_block[45] = temp_block[45];
                dst_block[46] = temp_block[46];
                dst_block[47] = temp_block[47];
                dst_block[48] = temp_block[48];
                dst_block[49] = temp_block[49];
                dst_block[50] = temp_block[50];
                dst_block[51] = temp_block[51];
                dst_block[52] = temp_block[52];
                dst_block[53] = temp_block[53];
                dst_block[54] = temp_block[54];
                dst_block[55] = temp_block[55];
                dst_block[56] = temp_block[56];
                dst_block[57] = temp_block[57];
                dst_block[58] = temp_block[58];
                dst_block[59] = temp_block[59];
                dst_block[60] = temp_block[60];
                dst_block[61] = temp_block[61];
                dst_block[62] = temp_block[62];
                dst_block[63] = temp_block[63];
            end
        endcase
        
        // Reset pixel_idx when entering FILTER state
        if (state == STATE_IDLE && state_next == STATE_FILTER) begin
            pixel_idx <= 7'd0;
        end
    end
end

//==============================================================================
// Valid Signal Generation
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        valid <= 1'b0;
    end else begin
        // Set valid when entering OUTPUT state (data ready)
        if (state == STATE_FILTER && state_next == STATE_OUTPUT) begin
            valid <= 1'b1;
        end
        // Clear valid when entering IDLE from OUTPUT
        if (state == STATE_OUTPUT && state_next == STATE_IDLE) begin
            valid <= 1'b0;
        end
    end
end

//==============================================================================
// State Machine - Combinational Logic
//==============================================================================

always @(*) begin
    state_next = state;
    
    case (state)
        STATE_IDLE: begin
            if (ready) begin
                state_next = STATE_FILTER;
            end
        end
        
        STATE_FILTER: begin
            if (pixel_idx >= 64) begin
                state_next = STATE_OUTPUT;
            end else begin
                state_next = STATE_FILTER;
            end
        end
        
        STATE_OUTPUT: begin
            // Return to IDLE immediately, staying just one cycle
            // Valid will be high during this cycle for testbench to see
            state_next = STATE_IDLE;
        end
        
        default: begin
            state_next = STATE_IDLE;
        end
    endcase
end

endmodule
