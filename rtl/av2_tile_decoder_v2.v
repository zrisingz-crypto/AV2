//==============================================================================
// AV2 Tile Decoder V2 - Uses fixed coefficient decoder
//==============================================================================

`timescale 1ns / 1ps

module av2_tile_decoder_v2 #(
    parameter MAX_WIDTH   = 128,
    parameter MAX_HEIGHT  = 128,
    parameter PIXEL_WIDTH = 10,
    parameter MAX_SB_SIZE = 64
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     start,
    
    input  wire [15:0]              frame_width,
    input  wire [15:0]              frame_height,
    input  wire [7:0]               qindex,
    input  wire [1:0]               frame_type,
    
    input  wire [127:0]             tile_data,
    input  wire                     tile_valid,
    output wire                     tile_ready,
    
    output reg  [31:0]              ref_read_addr,
    input  wire [9:0]               ref_pixel_data,
    output reg                      ref_read_en,
    
    output reg  [127:0]             recon_data,
    output reg  [31:0]              recon_addr,
    output reg                      recon_wr_en,
    output reg                      tile_done
);

// State machine
localparam IDLE             = 4'd0;
localparam PARSE_SB_HEADER  = 4'd1;
localparam ENTROPY_DECODE   = 4'd2;
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

// Block parameters
reg [5:0] block_width, block_height;
reg [5:0] block_x_coord, block_y_coord;
reg [6:0] intra_mode;
reg [3:0] tx_type;

// Loop variables
integer k, i, j;
reg [10:0] temp_recon;
reg [15:0] frame_idx;
reg [11:0] recon_pixel_idx;
reg [15:0] pixel_idx;

// Debug
reg [31:0] debug_cycle;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n)
        debug_cycle <= 0;
    else
        debug_cycle <= debug_cycle + 1;
end

// Entropy decoder
wire [15:0] entropy_symbol;
wire        entropy_valid;
reg         entropy_ready;
wire        entropy_done;

av2_entropy_decoder_real u_entropy (
    .clk(clk),
    .rst_n(rst_n),
    .bitstream_data(tile_data),
    .bitstream_valid(tile_valid && (state == ENTROPY_DECODE)),
    .bitstream_ready(tile_ready),
    .context_idx(16'd0),
    .context_prob(16'd16384),
    .context_idx_out(),
    .context_update_en(),
    .context_update_idx(),
    .context_update_bit(),
    .reset_contexts(state == IDLE && start),
    .symbol(entropy_symbol),
    .symbol_valid(entropy_valid),
    .symbol_ready(entropy_ready),
    .start(state == ENTROPY_DECODE),
    .done(entropy_done)
);

// Coefficient decoder (fixed version)
wire signed [15:0] coeff_out;
wire [11:0] coeff_addr;
wire coeff_valid;
reg coeff_ready;
wire [15:0] num_coeffs;
wire coeffs_valid;
reg coeffs_ready;
wire coeffs_done;

reg signed [15:0] decoded_coeffs[0:4095];

av2_coeff_decoder_fixed u_coeff (
    .clk(clk),
    .rst_n(rst_n),
    .context_idx(16'd0),
    .context_prob(16'd16384),
    .decoded_symbol(entropy_symbol),
    .symbol_valid(entropy_valid),
    .symbol_ready(entropy_ready),
    .coeff_out(coeff_out),
    .coeff_addr(coeff_addr),
    .coeff_valid(coeff_valid),
    .coeff_ready(coeff_ready),
    .num_coeffs(num_coeffs),
    .coeffs_valid(coeffs_valid),
    .coeffs_ready(coeffs_ready),
    .tx_size(block_width),
    .tx_type(tx_type),
    .qindex(qindex),
    .start(state == ENTROPY_DECODE),
    .done(coeffs_done)
);

// Store coefficients
always @(posedge clk) begin
    if (coeff_valid) begin
        coeff_ready <= 1'b1;
        if (coeff_addr < 4096)
            decoded_coeffs[coeff_addr] <= coeff_out;
    end else begin
        coeff_ready <= 1'b0;
    end
end

// Inverse transform
wire signed [15:0] residual_pixels[0:4095];
wire itx_valid;
reg itx_ready;
wire itx_done;

av2_inverse_transform_real_fixed u_itx (
    .clk(clk),
    .rst_n(rst_n),
    .coeffs(decoded_coeffs),
    .num_coeffs(num_coeffs),
    .tx_width(block_width),
    .tx_height(block_height),
    .tx_type(tx_type),
    .start(state == INVERSE_TX),
    .pixels(residual_pixels),
    .valid(itx_valid),
    .ready(itx_ready),
    .done(itx_done)
);

// Intra prediction
wire [9:0] intra_pred_pixels[0:4095];
wire intra_pred_valid;
reg intra_pred_ready;

reg [9:0] ref_top[0:63];
reg [9:0] ref_left[0:63];
reg [9:0] ref_top_left;

av2_intra_prediction_real_fixed u_pred (
    .clk(clk),
    .rst_n(rst_n),
    .ref_top(ref_top),
    .ref_left(ref_left),
    .ref_top_left(ref_top_left),
    .intra_mode(intra_mode),
    .block_width(block_width),
    .block_height(block_height),
    .start(state == PREDICTION && frame_type == 2'd0),
    .pred_pixels(intra_pred_pixels),
    .valid(intra_pred_valid),
    .ready(intra_pred_ready)
);

// Reconstruction buffers
reg [9:0] recon_frame[0:MAX_WIDTH*MAX_HEIGHT-1];
reg [9:0] pred_buffer[0:4095];

// Reference pixel update
always @(posedge clk) begin
    if (state == PARSE_SB_HEADER) begin
        if (block_y_coord > 0) begin
            for (k = 0; k < block_width; k = k + 1)
                ref_top[k] <= recon_frame[(block_y_coord - 1) * frame_width + block_x_coord + k];
        end else begin
            for (k = 0; k < block_width; k = k + 1)
                ref_top[k] <= 10'd128;
        end
        
        if (block_x_coord > 0) begin
            for (k = 0; k < block_height; k = k + 1)
                ref_left[k] <= recon_frame[(block_y_coord + k) * frame_width + block_x_coord - 1];
        end else begin
            for (k = 0; k < block_height; k = k + 1)
                ref_left[k] <= 10'd128;
        end
        
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
        coeffs_ready <= 0;
        itx_ready <= 0;
        intra_pred_ready <= 0;
    end else begin
        case (state)
            IDLE: begin
                tile_done <= 0;
                recon_wr_en <= 0;
                write_offset <= 0;
                coeffs_ready <= 0;
                itx_ready <= 0;
                intra_pred_ready <= 0;
                
                if (start) begin
                    sb_rows <= (frame_height + MAX_SB_SIZE - 1) / MAX_SB_SIZE;
                    sb_cols <= (frame_width + MAX_SB_SIZE - 1) / MAX_SB_SIZE;
                    sb_row <= 0;
                    sb_col <= 0;
                    total_pixels <= frame_width * frame_height;
                    
                    $display("[TIME %0t] V2: IDLE -> PARSE_SB_HEADER", $time);
                    state <= PARSE_SB_HEADER;
                end
            end
            
            PARSE_SB_HEADER: begin
                block_width <= 6'd16;
                block_height <= 6'd16;
                block_x_coord <= sb_col[5:0] * MAX_SB_SIZE;
                block_y_coord <= sb_row[5:0] * MAX_SB_SIZE;
                intra_mode <= 7'd0;
                tx_type <= 4'd0;
                
                if (debug_cycle < 100)
                    $display("[TIME %0t] V2: PARSE_SB_HEADER", $time);
                state <= ENTROPY_DECODE;
            end
            
            ENTROPY_DECODE: begin
                coeffs_ready <= 1'b1;
                
                if (debug_cycle < 100)
                    $display("[TIME %0t] V2: ENTROPY_DECODE, coeffs_done=%0b", $time, coeffs_done);
                
                if (coeffs_done) begin
                    coeffs_ready <= 1'b0;
                    $display("[TIME %0t] V2: ENTROPY_DECODE complete", $time);
                    state <= INVERSE_TX;
                end
            end
            
            INVERSE_TX: begin
                itx_ready <= 1'b1;
                
                if (debug_cycle < 100)
                    $display("[TIME %0t] V2: INVERSE_TX, itx_done=%0b", $time, itx_done);
                
                if (itx_done) begin
                    itx_ready <= 1'b0;
                    $display("[TIME %0t] V2: INVERSE_TX complete", $time);
                    state <= PREDICTION;
                end
            end
            
            PREDICTION: begin
                // Always set ready at the start of PREDICTION state
                intra_pred_ready <= 1'b1;
                
                if (debug_cycle < 200)
                    $display("[TIME %0t] V2: PREDICTION, intra_pred_valid=%0b", $time, intra_pred_valid);
                
                if (frame_type == 2'd0) begin
                    if (intra_pred_valid) begin
                        intra_pred_ready <= 1'b0;
                        $display("[TIME %0t] V2: PREDICTION complete", $time);
                        state <= RECONSTRUCTION;
                    end
                end else begin
                    intra_pred_ready <= 1'b0;
                    state <= RECONSTRUCTION;
                end
            end
            
            RECONSTRUCTION: begin
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
                
                for (j = 0; j < 16; j = j + 1) begin
                    recon_pixel_idx = j;
                    if (recon_pixel_idx < block_width * block_height) begin
                        frame_idx = (block_y_coord + (recon_pixel_idx / block_width)) * frame_width + 
                                    block_x_coord + (recon_pixel_idx % block_width);
                        if (frame_idx < MAX_WIDTH * MAX_HEIGHT)
                            recon_frame[frame_idx] <= pred_buffer[recon_pixel_idx];
                    end
                end
                
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
                    state <= DONE;
                end
            end
            
            DONE: begin
                recon_wr_en <= 0;
                tile_done <= 1;
                $display("[TIME %0t] V2: DONE!", $time);
                state <= IDLE;
            end
        endcase
    end
end

endmodule
