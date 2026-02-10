//==============================================================================
// AV2 Simplified Tile Decoder - For testing the pipeline without complex entropy decode
//==============================================================================

`timescale 1ns / 1ps

module av2_tile_decoder_simplified #(
    parameter MAX_WIDTH   = 128,
    parameter MAX_HEIGHT  = 128,
    parameter PIXEL_WIDTH = 10,
    parameter MAX_SB_SIZE = 64
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     start,
    
    // Frame parameters
    input  wire [15:0]              frame_width,
    input  wire [15:0]              frame_height,
    input  wire [7:0]               qindex,
    input  wire [1:0]               frame_type,
    
    // Bitstream input (not used in simplified version)
    input  wire [127:0]             tile_data,
    input  wire                     tile_valid,
    output wire                     tile_ready,
    
    // Reference frame interface
    output reg  [31:0]              ref_read_addr,
    input  wire [9:0]               ref_pixel_data,
    output reg                      ref_read_en,
    
    // Reconstructed frame output
    output reg  [127:0]             recon_data,
    output reg  [31:0]              recon_addr,
    output reg                      recon_wr_en,
    output reg                      tile_done
);

// State machine
localparam IDLE             = 4'd0;
localparam PARSE_SB_HEADER  = 4'd1;
localparam GEN_COEFFS       = 4'd2;  // Generate test coefficients
localparam INVERSE_TX       = 4'd3;
localparam PREDICTION       = 4'd4;
localparam RECONSTRUCTION   = 4'd5;
localparam CHECK_SB_COMPLETE = 4'd8;
localparam WRITE_OUTPUT     = 4'd6;
localparam DONE             = 4'd7;

reg [3:0] state;

// Superblock traversal
reg [15:0] sb_row, sb_col;
reg [15:0] sb_rows, sb_cols;
reg [15:0] write_offset;
reg [15:0] total_pixels;

// Current block parameters
reg [5:0] block_width, block_height;
reg [5:0] block_x_coord, block_y_coord;
reg [6:0] intra_mode;
reg [3:0] tx_type;

// Counters for sub-modules
reg [7:0] itx_counter;
reg [7:0] pred_counter;

// Loop variables
integer k, i, j;

// Temporary variables
reg [10:0] temp_recon;
reg [15:0] frame_idx;
reg [11:0] recon_pixel_idx;
reg [15:0] pixel_idx;

// Debug cycle counter
reg [31:0] debug_cycle;

assign tile_ready = 1'b1;  // Always ready in simplified version

// Test coefficients and residual pixels
reg signed [15:0] decoded_coeffs[0:4095];
reg signed [15:0] residual_pixels[0:4095];

// Predicted pixels
reg [9:0]  intra_pred_pixels[0:4095];

// Reference pixels
reg [9:0]  ref_top[0:63];
reg [9:0]  ref_left[0:63];
reg [9:0]  ref_top_left;

// Reconstruction buffer
reg [9:0] recon_frame[0:MAX_WIDTH*MAX_HEIGHT-1];
reg [9:0] pred_buffer[0:4095];

// Debug cycle counter
always @(posedge clk or negedge rst_n) begin
    if (!rst_n)
        debug_cycle <= 32'd0;
    else
        debug_cycle <= debug_cycle + 1;
end

// Reference pixel update
always @(posedge clk) begin
    if (state == PARSE_SB_HEADER) begin
        // Update top reference pixels
        if (block_y_coord > 0) begin
            for (k = 0; k < block_width; k = k + 1)
                ref_top[k] <= recon_frame[(block_y_coord - 1) * frame_width + block_x_coord + k];
        end else begin
            for (k = 0; k < block_width; k = k + 1)
                ref_top[k] <= 10'd128;
        end
        
        // Update left reference pixels
        if (block_x_coord > 0) begin
            for (k = 0; k < block_height; k = k + 1)
                ref_left[k] <= recon_frame[(block_y_coord + k) * frame_width + block_x_coord - 1];
        end else begin
            for (k = 0; k < block_height; k = k + 1)
                ref_left[k] <= 10'd128;
        end
        
        // Update top-left reference pixel
        if (block_x_coord > 0 && block_y_coord > 0)
            ref_top_left <= recon_frame[(block_y_coord - 1) * frame_width + block_x_coord - 1];
        else
            ref_top_left <= 10'd128;
    end
end

// Main state machine
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        sb_row <= 0;
        sb_col <= 0;
        tile_done <= 0;
        recon_wr_en <= 0;
        write_offset <= 0;
        itx_counter <= 0;
        pred_counter <= 0;
    end else begin
        case (state)
            IDLE: begin
                tile_done <= 0;
                recon_wr_en <= 0;
                write_offset <= 0;
                itx_counter <= 0;
                pred_counter <= 0;
                
                if (start) begin
                    sb_rows <= (frame_height + MAX_SB_SIZE - 1) / MAX_SB_SIZE;
                    sb_cols <= (frame_width + MAX_SB_SIZE - 1) / MAX_SB_SIZE;
                    sb_row <= 0;
                    sb_col <= 0;
                    total_pixels <= frame_width * frame_height;
                    
                    $display("[TIME %0t] IDLE -> PARSE_SB_HEADER, frame=%0dx%0d", 
                             $time, frame_width, frame_height);
                    state <= PARSE_SB_HEADER;
                end
            end
            
            PARSE_SB_HEADER: begin
                block_width <= 6'd16;
                block_height <= 6'd16;
                block_x_coord <= sb_col[5:0] * MAX_SB_SIZE;
                block_y_coord <= sb_row[5:0] * MAX_SB_SIZE;
                intra_mode <= 7'd0;  // DC mode
                tx_type <= 4'd0;     // DCT_DCT
                
                $display("[TIME %0t] PARSE_SB_HEADER: sb_row=%0d, sb_col=%0d", 
                         $time, sb_row, sb_col);
                state <= GEN_COEFFS;
            end
            
            GEN_COEFFS: begin
                // Generate test coefficients (simplified)
                for (k = 0; k < 256; k = k + 1) begin
                    decoded_coeffs[k] <= 16'sd10;  // Small residual
                end
                
                // Generate residual pixels directly (bypass inverse transform)
                for (k = 0; k < 256; k = k + 1) begin
                    residual_pixels[k] <= 16'sd10;
                end
                
                $display("[TIME %0t] GEN_COEFFS -> INVERSE_TX", $time);
                state <= INVERSE_TX;
            end
            
            INVERSE_TX: begin
                // Simulate inverse transform delay
                itx_counter <= itx_counter + 1;
                if (itx_counter >= 2) begin
                    itx_counter <= 0;
                    $display("[TIME %0t] INVERSE_TX -> PREDICTION", $time);
                    state <= PREDICTION;
                end
            end
            
            PREDICTION: begin
                // Generate prediction (DC mode - simplified)
                pred_counter <= pred_counter + 1;
                if (pred_counter >= 2) begin
                    pred_counter <= 0;
                    
                    // Simple DC prediction
                    for (k = 0; k < 256; k = k + 1) begin
                        intra_pred_pixels[k] <= 10'd128;
                    end
                    
                    $display("[TIME %0t] PREDICTION -> RECONSTRUCTION", $time);
                    state <= RECONSTRUCTION;
                end
            end
            
            RECONSTRUCTION: begin
                // Add prediction and residual with clipping
                for (i = 0; i < 16; i = i + 1) begin
                    recon_pixel_idx = i;
                    if (recon_pixel_idx < block_width * block_height) begin
                        temp_recon = intra_pred_pixels[recon_pixel_idx] + 
                                     residual_pixels[recon_pixel_idx];
                        
                        if (temp_recon > 10'd1023)
                            pred_buffer[recon_pixel_idx] <= 10'd1023;
                        else if (temp_recon < 10'd0)
                            pred_buffer[recon_pixel_idx] <= 10'd0;
                        else
                            pred_buffer[recon_pixel_idx] <= temp_recon[9:0];
                    end
                end
                
                // Store in reconstruction frame
                for (j = 0; j < 16; j = j + 1) begin
                    recon_pixel_idx = j;
                    if (recon_pixel_idx < block_width * block_height) begin
                        frame_idx = (block_y_coord + (recon_pixel_idx / block_width)) * frame_width + 
                                    block_x_coord + (recon_pixel_idx % block_width);
                        if (frame_idx < MAX_WIDTH * MAX_HEIGHT)
                            recon_frame[frame_idx] <= pred_buffer[recon_pixel_idx];
                    end
                end
                
                $display("[TIME %0t] RECONSTRUCTION -> CHECK_SB_COMPLETE", $time);
                state <= CHECK_SB_COMPLETE;
            end
            
            CHECK_SB_COMPLETE: begin
                if (sb_col < sb_cols - 1) begin
                    sb_col <= sb_col + 1;
                    state <= PARSE_SB_HEADER;
                end else if (sb_row < sb_rows - 1) begin
                    sb_col <= 0;
                    sb_row <= sb_row + 1;
                    state <= PARSE_SB_HEADER;
                end else begin
                    $display("[TIME %0t] CHECK_SB_COMPLETE -> WRITE_OUTPUT", $time);
                    state <= WRITE_OUTPUT;
                end
            end
            
            WRITE_OUTPUT: begin
                recon_wr_en <= 1;
                
                if (write_offset < total_pixels) begin
                    recon_addr <= write_offset;
                    
                    for (k = 0; k < 16; k = k + 1) begin
                        pixel_idx = write_offset + k;
                        case (k)
                            4'd0:  recon_data[7:0]   <= recon_frame[pixel_idx];
                            4'd1:  recon_data[15:8]  <= recon_frame[pixel_idx];
                            4'd2:  recon_data[23:16] <= recon_frame[pixel_idx];
                            4'd3:  recon_data[31:24] <= recon_frame[pixel_idx];
                            4'd4:  recon_data[39:32] <= recon_frame[pixel_idx];
                            4'd5:  recon_data[47:40] <= recon_frame[pixel_idx];
                            4'd6:  recon_data[55:48] <= recon_frame[pixel_idx];
                            4'd7:  recon_data[63:56] <= recon_frame[pixel_idx];
                            4'd8:  recon_data[71:64] <= recon_frame[pixel_idx];
                            4'd9:  recon_data[79:72] <= recon_frame[pixel_idx];
                            4'd10: recon_data[87:80] <= recon_frame[pixel_idx];
                            4'd11: recon_data[95:88] <= recon_frame[pixel_idx];
                            4'd12: recon_data[103:96] <= recon_frame[pixel_idx];
                            4'd13: recon_data[111:104] <= recon_frame[pixel_idx];
                            4'd14: recon_data[119:112] <= recon_frame[pixel_idx];
                            4'd15: recon_data[127:120] <= recon_frame[pixel_idx];
                        endcase
                    end
                    
                    write_offset <= write_offset + 16;
                end else begin
                    $display("[TIME %0t] WRITE_OUTPUT -> DONE", $time);
                    state <= DONE;
                end
            end
            
            DONE: begin
                recon_wr_en <= 0;
                tile_done <= 1;
                $display("[TIME %0t] DONE!", $time);
                state <= IDLE;
            end
        endcase
    end
end

endmodule
