//==============================================================================
// AV2 Real Tile Decoder - Integrated with Real Decode Modules (FIXED VERSION)
// Implements full decode pipeline with real: Entropy Decode → Inverse TX → Prediction → Reconstruction
//==============================================================================

`timescale 1ns / 1ps

module av2_tile_decoder_real_fixed #(
    parameter MAX_WIDTH   = 128,
    parameter MAX_HEIGHT  = 128,
    parameter PIXEL_WIDTH = 10,
    parameter MAX_SB_SIZE = 64  // Superblock size
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     start,
    
    // Frame parameters
    input  wire [15:0]              frame_width,
    input  wire [15:0]              frame_height,
    input  wire [7:0]               qindex,
    input  wire [1:0]               frame_type,  // 0:KEY, 1:INTER
    
    // Bitstream input
    input  wire [127:0]             tile_data,
    input  wire                     tile_valid,
    output wire                     tile_ready,
    
    // Reference frame interface (for inter prediction)
    output reg  [31:0]              ref_read_addr,
    input  wire [9:0]               ref_pixel_data,
    output reg                      ref_read_en,
    
    // Reconstructed frame output
    output reg  [127:0]             recon_data,
    output reg  [31:0]              recon_addr,
    output reg                      recon_wr_en,
    output reg                      tile_done
);

//==============================================================================
// Decoding State Machine
//==============================================================================

localparam IDLE             = 4'd0;
localparam PARSE_SB_HEADER  = 4'd1;  // Parse superblock header
localparam ENTROPY_DECODE   = 4'd2;  // Entropy decode
localparam INVERSE_TX       = 4'd3;  // Inverse transform
localparam PREDICTION       = 4'd4;  // Prediction
localparam RECONSTRUCTION   = 4'd5;  // Reconstruction
localparam CHECK_SB_COMPLETE = 4'd8; // Check superblock completion
localparam WRITE_OUTPUT     = 4'd6;  // Write output
localparam DONE             = 4'd7;

reg [3:0] state, state_next;

// Superblock traversal
reg [15:0] sb_row, sb_col;
reg [15:0] sb_rows, sb_cols;

// Bitstream consumption tracking
reg [31:0] bitstream_count;
reg        bitstream_complete;

// Current block parameters
reg [5:0] block_width, block_height;
reg [5:0] block_x_coord, block_y_coord;  // Block coordinates
reg [6:0] intra_mode;
reg signed [15:0] mv_x, mv_y;
reg [3:0] tx_type;

// Output write control
reg [15:0] write_offset;
reg [15:0] total_pixels;

// Loop variables - must be declared at module level
integer k, i, j;

// Temporary variables for RECONSTRUCTION state (declared at module level)
reg [10:0] temp_recon;
reg [15:0] frame_idx;
reg [11:0] recon_pixel_idx;

// Temporary variables for WRITE_OUTPUT state
reg [15:0] pixel_idx;

// Debug cycle counter
reg [31:0] debug_cycle;

//==============================================================================
// Real Module Instantiations
//==============================================================================

// 1. Real Entropy Decoder
wire [15:0] entropy_symbol;
wire        entropy_valid;
wire        entropy_ready;
wire        entropy_done;

wire [15:0] context_idx_out;
wire        context_update_en;
wire [15:0] context_update_idx;
wire        context_update_bit;

av2_entropy_decoder_real #(
    .DATA_WIDTH(128)
) u_entropy_decoder_real (
    .clk              (clk),
    .rst_n            (rst_n),
    .bitstream_data   (tile_data),
    .bitstream_valid  (tile_valid && (state == ENTROPY_DECODE)),
    .bitstream_ready  (tile_ready),
    .context_idx      (16'd0),  // Simplified: fixed context
    .context_prob     (16'd16384),  // 50% probability
    .context_idx_out  (context_idx_out),
    .context_update_en (context_update_en),
    .context_update_idx(context_update_idx),
    .context_update_bit(context_update_bit),
    .reset_contexts   (state == IDLE && start),
    .symbol           (entropy_symbol),
    .symbol_valid     (entropy_valid),
    .symbol_ready     (entropy_ready),
    .start            (state == ENTROPY_DECODE),
    .done             (entropy_done)
);

// Local control for entropy decoder
reg local_entropy_ready;
assign entropy_ready = local_entropy_ready;

// 2. Real Coefficient Decoder
reg signed [15:0] decoded_coeffs[0:4095];  // Max 64x64 block = 4096 coefficients
wire [15:0]        num_coeffs_out;
wire               coeffs_valid;
reg                coeffs_ready;
wire               coeffs_done;

// For individual coefficient output
reg signed [15:0]  coeff_out_temp;
reg [11:0]         coeff_addr_temp;
wire               coeff_valid_temp;
reg                coeff_ready_temp;

av2_coeff_decoder_real #(
    .MAX_COEFFS(4096),
    .MAX_TX_SIZE(64)
) u_coeff_decoder_real (
    .clk              (clk),
    .rst_n            (rst_n),
    .context_idx      (context_idx_out),
    .context_prob     (16'd16384),  // Simplified: fixed probability
    .decoded_symbol   (entropy_symbol),
    .symbol_valid     (entropy_valid),
    .symbol_ready     (entropy_ready),
    .coeff_out        (coeff_out_temp),
    .coeff_addr       (coeff_addr_temp),
    .coeff_valid      (coeff_valid_temp),
    .coeff_ready      (coeff_ready_temp),
    .num_coeffs       (num_coeffs_out),
    .coeffs_valid     (coeffs_valid),
    .coeffs_ready     (coeffs_ready),
    .tx_size          (block_width),
    .tx_type          (tx_type),
    .qindex           (qindex),
    .start            (state == ENTROPY_DECODE),
    .done             (coeffs_done)
);

// Transfer coefficients from coeff decoder to internal memory
always @(posedge clk) begin
    if (coeffs_valid && coeff_valid_temp) begin
        coeff_ready_temp <= 1'b1;
        if (coeff_addr_temp < num_coeffs_out) begin
            decoded_coeffs[coeff_addr_temp] <= coeff_out_temp;
        end
    end else begin
        coeff_ready_temp <= 1'b0;
    end
end

// 3. Real Inverse Transform
wire signed [15:0] residual_pixels[0:4095];
wire               itx_valid;
reg                itx_ready;
wire               itx_done;

av2_inverse_transform_real_fixed #(
    .MAX_TX_SIZE(64),
    .BIT_DEPTH(10)
) u_inverse_transform_real (
    .clk              (clk),
    .rst_n            (rst_n),
    .coeffs           (decoded_coeffs),
    .num_coeffs       (num_coeffs_out),
    .tx_width         (block_width),
    .tx_height        (block_height),
    .tx_type          (tx_type),
    .start            (state == INVERSE_TX),
    .pixels           (residual_pixels),
    .valid            (itx_valid),
    .ready            (itx_ready),
    .done             (itx_done)
);

// 4. Real Intra Prediction
wire [9:0]  intra_pred_pixels[0:4095];
wire        intra_pred_valid;
reg         intra_pred_ready;
wire        intra_pred_done;

// Reference pixels
reg [9:0]  ref_top[0:63];
reg [9:0]  ref_left[0:63];
reg [9:0]  ref_top_left;

av2_intra_prediction_real_fixed #(
    .MAX_BLOCK_SIZE(64),
    .BIT_DEPTH(10)
) u_intra_prediction_real (
    .clk              (clk),
    .rst_n            (rst_n),
    .ref_top          (ref_top),
    .ref_left         (ref_left),
    .ref_top_left     (ref_top_left),
    .intra_mode       (intra_mode),
    .block_width      (block_width),
    .block_height     (block_height),
    .start            (state == PREDICTION && frame_type == 2'd0),
    .pred_pixels      (intra_pred_pixels),
    .valid            (intra_pred_valid),
    .ready            (intra_pred_ready)
);

//==============================================================================
// Reconstruction Buffer
//==============================================================================

reg [9:0] recon_frame[0:MAX_WIDTH*MAX_HEIGHT-1];
reg [9:0] pred_buffer[0:4095];

//==============================================================================
// Reference Pixel Management
//==============================================================================

// Update reference pixels
always @(posedge clk) begin
    if (state == PARSE_SB_HEADER) begin
        // Update top reference pixels
        if (block_y_coord > 0) begin
            for (k = 0; k < block_width; k = k + 1) begin
                ref_top[k] <= recon_frame[(block_y_coord - 1) * frame_width + block_x_coord + k];
            end
        end else begin
            // Boundary: use default values
            for (k = 0; k < block_width; k = k + 1) begin
                ref_top[k] <= 10'd128;
            end
        end
        
        // Update left reference pixels
        if (block_x_coord > 0) begin
            for (k = 0; k < block_height; k = k + 1) begin
                ref_left[k] <= recon_frame[(block_y_coord + k) * frame_width + block_x_coord - 1];
            end
        end else begin
            // Boundary: use default values
            for (k = 0; k < block_height; k = k + 1) begin
                ref_left[k] <= 10'd128;
            end
        end
        
        // Update top-left reference pixel
        if (block_x_coord > 0 && block_y_coord > 0) begin
            ref_top_left <= recon_frame[(block_y_coord - 1) * frame_width + block_x_coord - 1];
        end else begin
            ref_top_left <= 10'd128;
        end
    end
end

//==============================================================================
// Debug Cycle Counter
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        debug_cycle <= 32'd0;
    end else begin
        debug_cycle <= debug_cycle + 1;
    end
end

//==============================================================================
// Main Control Logic - Sequential
//==============================================================================

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        sb_row <= 16'd0;
        sb_col <= 16'd0;
        tile_done <= 1'b0;
        recon_wr_en <= 1'b0;
        itx_ready <= 1'b0;
        intra_pred_ready <= 1'b0;
        local_entropy_ready <= 1'b0;
        
        for (k = 0; k < 4096; k = k + 1) begin
            pred_buffer[k] <= 10'd0;
        end
    end else begin
        state <= state_next;
        
        case (state)
            IDLE: begin
                tile_done <= 1'b0;
                recon_wr_en <= 1'b0;
                write_offset <= 16'd0;
                bitstream_count <= 32'd0;
                bitstream_complete <= 1'b0;
                
                if (start) begin
                    // Calculate superblock count
                    sb_rows <= (frame_height + MAX_SB_SIZE - 1) / MAX_SB_SIZE;
                    sb_cols <= (frame_width + MAX_SB_SIZE - 1) / MAX_SB_SIZE;
                    sb_row <= 16'd0;
                    sb_col <= 16'd0;
                    
                    // Calculate total pixels for output
                    total_pixels <= frame_width * frame_height;
                    
                    if (debug_cycle < 100)
                        $display("[TIME %0t] IDLE -> START, frame=%0dx%0d, sb_rows=%0d, sb_cols=%0d", 
                                 $time, frame_width, frame_height, 
                                 (frame_height + MAX_SB_SIZE - 1) / MAX_SB_SIZE,
                                 (frame_width + MAX_SB_SIZE - 1) / MAX_SB_SIZE);
                end
            end
            
            PARSE_SB_HEADER: begin
                // Parse superblock header (simplified: fixed parameters)
                block_width <= 6'd16;
                block_height <= 6'd16;
                block_x_coord <= sb_col[5:0] * MAX_SB_SIZE;
                block_y_coord <= sb_row[5:0] * MAX_SB_SIZE;
                
                if (frame_type == 2'd0)
                    intra_mode <= 7'd0;  // DC mode
                else begin
                    mv_x <= 16'sd0;
                    mv_y <= 16'sd0;
                end
                
                tx_type <= 4'd0;  // DCT_DCT
                
                if (debug_cycle < 500)
                    $display("[TIME %0t] PARSE_SB_HEADER: sb_row=%0d, sb_col=%0d", 
                             $time, sb_row, sb_col);
            end
            
            ENTROPY_DECODE: begin
                // Wait for real entropy decoder and coefficient decoder
                local_entropy_ready <= 1'b1;
                coeffs_ready <= 1'b1;
                
                if (debug_cycle < 100)
                    $display("[TIME %0t] ENTROPY_DECODE: entropy_done=%0b, coeffs_done=%0b", 
                             $time, entropy_done, coeffs_done);
                
                if (coeffs_done) begin
                    local_entropy_ready <= 1'b0;
                    coeffs_ready <= 1'b0;
                    $display("[TIME %0t] ENTROPY_DECODE complete", $time);
                end
            end
            
            INVERSE_TX: begin
                // Wait for real inverse transform
                itx_ready <= 1'b1;
                if (itx_done) begin
                    itx_ready <= 1'b0;
                    
                    if (debug_cycle < 500)
                        $display("[TIME %0t] INVERSE_TX complete", $time);
                end
            end
            
            PREDICTION: begin
                // Wait for real intra prediction
                if (frame_type == 2'd0) begin
                    intra_pred_ready <= 1'b1;
                    if (intra_pred_valid) begin
                        intra_pred_ready <= 1'b0;
                        
                        if (debug_cycle < 500)
                            $display("[TIME %0t] PREDICTION complete", $time);
                    end
                end
            end
            
            RECONSTRUCTION: begin
                // Real reconstruction: add prediction and residual
                // Process one row at a time to avoid combinational loops
                for (i = 0; i < 16; i = i + 1) begin
                    recon_pixel_idx = i;
                    if (recon_pixel_idx < block_width * block_height) begin
                        // pred + residual
                        temp_recon = intra_pred_pixels[recon_pixel_idx] + residual_pixels[recon_pixel_idx];
                        
                        // Clipping
                        if (temp_recon > 10'd1023) begin
                            pred_buffer[recon_pixel_idx] <= 10'd1023;
                        end else if (temp_recon < 10'd0) begin
                            pred_buffer[recon_pixel_idx] <= 10'd0;
                        end else begin
                            pred_buffer[recon_pixel_idx] <= temp_recon[9:0];
                        end
                    end
                end
                
                // Store in reconstruction frame - process in subsequent cycles
                for (j = 0; j < 16; j = j + 1) begin
                    recon_pixel_idx = j;
                    if (recon_pixel_idx < block_width * block_height) begin
                        frame_idx = (block_y_coord + (recon_pixel_idx / block_width)) * frame_width + 
                                    block_x_coord + (recon_pixel_idx % block_width);
                        if (frame_idx < MAX_WIDTH * MAX_HEIGHT) begin
                            recon_frame[frame_idx] <= pred_buffer[recon_pixel_idx];
                        end
                    end
                end
                
                local_entropy_ready <= 1'b0;
                itx_ready <= 1'b0;
                intra_pred_ready <= 1'b0;
                
                if (debug_cycle < 500)
                    $display("[TIME %0t] RECONSTRUCTION done", $time);
            end
            
            CHECK_SB_COMPLETE: begin
                // Update superblock counters
                if (sb_col < sb_cols - 1) begin
                    // Move to next superblock in same row
                    sb_col <= sb_col + 16'd1;
                    
                    if (debug_cycle < 500)
                        $display("[TIME %0t] CHECK_SB_COMPLETE: next col, sb_col=%0d", $time, sb_col + 1);
                end else if (sb_row < sb_rows - 1) begin
                    // Move to next row, reset column
                    sb_col <= 16'd0;
                    sb_row <= sb_row + 16'd1;
                    
                    if (debug_cycle < 500)
                        $display("[TIME %0t] CHECK_SB_COMPLETE: next row, sb_row=%0d", $time, sb_row + 1);
                end else begin
                    if (debug_cycle < 500)
                        $display("[TIME %0t] CHECK_SB_COMPLETE: all SBs done, goto WRITE_OUTPUT", $time);
                end
            end
            
            WRITE_OUTPUT: begin
                // Write reconstructed data to frame buffer
                recon_wr_en <= 1'b1;
                
                if (write_offset < total_pixels) begin
                    recon_addr <= write_offset;
                    
                    // Pack 16 pixels into 128-bit word
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
                    
                    write_offset <= write_offset + 16'd16;
                    
                    if (debug_cycle < 200 && (write_offset % 256 == 0))
                        $display("[TIME %0t] WRITE_OUTPUT: offset=%0d / %0d", $time, write_offset, total_pixels);
                end else begin
                    if (debug_cycle < 500)
                        $display("[TIME %0t] WRITE_OUTPUT complete, offset=%0d >= %0d", $time, write_offset, total_pixels);
                end
            end
            
            DONE: begin
                recon_wr_en <= 1'b0;
                tile_done <= 1'b1;
                
                if (debug_cycle < 500)
                    $display("[TIME %0t] DONE asserted", $time);
            end
        endcase
    end
end

//==============================================================================
// State Transition Logic - Combinational
//==============================================================================
always @(*) begin
    state_next = state;
    
    case (state)
        IDLE: begin
            if (start)
                state_next = PARSE_SB_HEADER;
        end
        
        PARSE_SB_HEADER: begin
            state_next = ENTROPY_DECODE;
        end
        
        ENTROPY_DECODE: begin
            if (coeffs_done)
                state_next = INVERSE_TX;
        end
        
        INVERSE_TX: begin
            if (itx_done)
                state_next = PREDICTION;
        end
        
        PREDICTION: begin
            if (frame_type == 2'd0) begin
                if (intra_pred_valid)
                    state_next = RECONSTRUCTION;
            end else begin
                // Inter frame - skip prediction for now
                state_next = RECONSTRUCTION;
            end
        end
        
        RECONSTRUCTION: begin
            // Single cycle reconstruction, then check SB
            state_next = CHECK_SB_COMPLETE;
        end
        
        CHECK_SB_COMPLETE: begin
            // Check if current superblock row is complete
            if (sb_col < sb_cols - 1) begin
                // Move to next superblock in same row
                state_next = PARSE_SB_HEADER;
            end else if (sb_row < sb_rows - 1) begin
                // Move to next row
                state_next = PARSE_SB_HEADER;
            end else begin
                // All superblocks complete, write output
                state_next = WRITE_OUTPUT;
            end
        end
        
        WRITE_OUTPUT: begin
            // Use >= to ensure we catch the completion
            if (write_offset >= total_pixels)
                state_next = DONE;
        end
        
        DONE: begin
            state_next = IDLE;
        end
        
        default: begin
            state_next = IDLE;
        end
    endcase
end

endmodule
