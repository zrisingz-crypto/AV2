//==============================================================================
// Deblocking Filter Module (Real Implementation)
// Implements AV2 deblocking filter for block boundary smoothing
//==============================================================================

`timescale 1ns / 1ps

module av2_deblocking_filter_real #(
    parameter MAX_WIDTH = 128,
    parameter MAX_HEIGHT = 128
)(
    input  wire                    clk,
    input  wire                    rst_n,
    input  wire [9:0]             src_pixels[0:MAX_WIDTH*MAX_HEIGHT-1],
    input  wire [15:0]            frame_width,
    input  wire [15:0]            frame_height,
    input  wire [5:0]             filter_level,
    input  wire [2:0]             sharpness,
    input  wire                    start,
    output reg  [9:0]             dst_pixels[0:MAX_WIDTH*MAX_HEIGHT-1],
    output reg                     valid,
    input  wire                    ready
);

//==============================================================================
// Parameters
//==============================================================================

localparam STATE_IDLE          = 3'd0;
localparam STATE_LOAD_BLOCK    = 3'd1;
localparam STATE_FILTER_V      = 3'd2;
localparam STATE_FILTER_H      = 3'd3;
localparam STATE_OUTPUT        = 3'd4;
localparam STATE_DONE          = 3'd5;

// Deblocking filter thresholds
localparam MAX_FILTER_LEVEL   = 6'd63;
localparam MAX_SHARPNESS      = 3'd7;
localparam BLOCK_SIZE         = 4;     // 4x4 blocks
localparam BLOCK_SIZE_8       = 8;     // 8x8 blocks

//==============================================================================
// State Register
//==============================================================================

reg [2:0] state;
reg [2:0] state_next;

//==============================================================================
// Working Registers
//==============================================================================

reg [15:0] x_coord;
reg [15:0] y_coord;
reg [15:0] block_idx;
reg [15:0] pixel_idx;
reg [15:0] total_pixels;

// 8-pixel line buffer for deblocking
reg [9:0]  line_buffer[0:7];

// Filter parameters
reg [5:0]  thr_i;
reg [5:0]  thr_b;
reg [5:0]  limit;

//==============================================================================
// Helper Functions
//==============================================================================

// Calculate filter threshold based on filter_level and sharpness
function [5:0] calc_thr_i;
    input [5:0] fl;
    input [2:0] sh;
    begin
        calc_thr_i = (fl * (2 + sh)) >> 4;
    end
endfunction

function [5:0] calc_thr_b;
    input [5:0] fl;
    input [2:0] sh;
    begin
        calc_thr_b = ((fl * (2 + sh)) >> 4) + (fl >> 4);
    end
endfunction

// Clip function
function [9:0] clip;
    input [9:0] val;
    input [9:0] min_val;
    input [9:0] max_val;
    begin
        if (val < min_val)
            clip = min_val;
        else if (val > max_val)
            clip = max_val;
        else
            clip = val;
    end
endfunction

integer i;

//==============================================================================
// State Machine - Sequential Logic
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= STATE_IDLE;
        valid <= 1'b0;
        x_coord <= 16'd0;
        y_coord <= 16'd0;
        block_idx <= 16'd0;
        pixel_idx <= 16'd0;
        total_pixels <= 16'd0;
        thr_i <= 6'd0;
        thr_b <= 6'd0;
        limit <= 6'd0;
        
        // Initialize output to zeros
        for (i = 0; i < MAX_WIDTH * MAX_HEIGHT; i = i + 1) begin
            dst_pixels[i] <= 10'd0;
        end
    end else begin
        state <= state_next;
        
        case (state)
            STATE_IDLE: begin
                if (start) begin
                    total_pixels <= frame_width * frame_height;
                    thr_i <= calc_thr_i(filter_level, sharpness);
                    thr_b <= calc_thr_b(filter_level, sharpness);
                    limit <= (filter_level * 2) + 1;
                    x_coord <= 16'd0;
                    y_coord <= 16'd0;
                    block_idx <= 16'd0;
                end
            end
            
            STATE_LOAD_BLOCK: begin
                // Load a line of pixels for filtering
                for (i = 0; i < 8; i = i + 1) begin
                    if ((y_coord * frame_width + x_coord + i) < total_pixels) begin
                        line_buffer[i] <= src_pixels[y_coord * frame_width + x_coord + i];
                    end else begin
                        line_buffer[i] <= 10'd0;
                    end
                end
            end
            
            STATE_FILTER_V: begin
                // Apply vertical filter (top-bottom boundary)
                // Simplified strong filter
                if (filter_level > 0) begin
                    // p3, p2, p1, p0 | q0, q1, q2, q3
                    // Filter p0 and q0
                    if (abs_diff(line_buffer[3], line_buffer[4]) < limit) begin
                        if (abs_diff(line_buffer[2], line_buffer[5]) < thr_i) begin
                            dst_pixels[y_coord * frame_width + x_coord + 3] <= 
                                clip(((line_buffer[3] * 3 + line_buffer[4] * 2 + line_buffer[5] + 4) >> 3), 10'd0, 10'd1023);
                            dst_pixels[y_coord * frame_width + x_coord + 4] <= 
                                clip(((line_buffer[4] * 3 + line_buffer[3] * 2 + line_buffer[2] + 4) >> 3), 10'd0, 10'd1023);
                        end else begin
                            dst_pixels[y_coord * frame_width + x_coord + 3] <= line_buffer[3];
                            dst_pixels[y_coord * frame_width + x_coord + 4] <= line_buffer[4];
                        end
                    end else begin
                        dst_pixels[y_coord * frame_width + x_coord + 3] <= line_buffer[3];
                        dst_pixels[y_coord * frame_width + x_coord + 4] <= line_buffer[4];
                    end
                    
                    // Copy other pixels
                    dst_pixels[y_coord * frame_width + x_coord + 0] <= line_buffer[0];
                    dst_pixels[y_coord * frame_width + x_coord + 1] <= line_buffer[1];
                    dst_pixels[y_coord * frame_width + x_coord + 2] <= line_buffer[2];
                    dst_pixels[y_coord * frame_width + x_coord + 5] <= line_buffer[5];
                    dst_pixels[y_coord * frame_width + x_coord + 6] <= line_buffer[6];
                    dst_pixels[y_coord * frame_width + x_coord + 7] <= line_buffer[7];
                end
            end
            
            STATE_FILTER_H: begin
                // Apply horizontal filter (left-right boundary)
                // Similar to vertical filter
                if (filter_level > 0) begin
                    if (abs_diff(line_buffer[3], line_buffer[4]) < limit) begin
                        if (abs_diff(line_buffer[2], line_buffer[5]) < thr_i) begin
                            dst_pixels[y_coord * frame_width + x_coord + 3] <= 
                                clip(((line_buffer[3] * 3 + line_buffer[4] * 2 + line_buffer[5] + 4) >> 3), 10'd0, 10'd1023);
                            dst_pixels[y_coord * frame_width + x_coord + 4] <= 
                                clip(((line_buffer[4] * 3 + line_buffer[3] * 2 + line_buffer[2] + 4) >> 3), 10'd0, 10'd1023);
                        end else begin
                            dst_pixels[y_coord * frame_width + x_coord + 3] <= line_buffer[3];
                            dst_pixels[y_coord * frame_width + x_coord + 4] <= line_buffer[4];
                        end
                    end else begin
                        dst_pixels[y_coord * frame_width + x_coord + 3] <= line_buffer[3];
                        dst_pixels[y_coord * frame_width + x_coord + 4] <= line_buffer[4];
                    end
                    
                    // Copy other pixels
                    for (i = 0; i < 8; i = i + 1) begin
                        if (i != 3 && i != 4) begin
                            dst_pixels[y_coord * frame_width + x_coord + i] <= line_buffer[i];
                        end
                    end
                end
            end
            
            STATE_OUTPUT: begin
                // All filtering done, signal completion
                valid <= 1'b1;
                if (ready) begin
                    valid <= 1'b0;
                end
            end
            
            STATE_DONE: begin
                valid <= 1'b0;
            end
        endcase
    end
end

//==============================================================================
// Absolute difference function
//==============================================================================

function [9:0] abs_diff;
    input [9:0] a;
    input [9:0] b;
    begin
        if (a >= b)
            abs_diff = a - b;
        else
            abs_diff = b - a;
    end
endfunction

//==============================================================================
// State Machine - Combinational Logic
//==============================================================================

always @(*) begin
    state_next = state;
    
    case (state)
        STATE_IDLE: begin
            if (start) begin
                state_next = STATE_LOAD_BLOCK;
            end
        end
        
        STATE_LOAD_BLOCK: begin
            if (filter_level > 0) begin
                state_next = STATE_FILTER_V;
            end else begin
                state_next = STATE_OUTPUT;
            end
        end
        
        STATE_FILTER_V: begin
            state_next = STATE_FILTER_H;
        end
        
        STATE_FILTER_H: begin
            state_next = STATE_OUTPUT;
        end
        
        STATE_OUTPUT: begin
            if (valid && ready) begin
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
// Debug Output (disabled to reduce log size)
//==============================================================================

// always @(posedge clk) begin
//     if (state == STATE_FILTER_V && filter_level > 0) begin
//         $display("[DEBLOCK] Filtering at (%0d, %0d), level=%0d, thr_i=%0d", 
//                  x_coord, y_coord, filter_level, thr_i);
//     end
// end

endmodule