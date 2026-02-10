//==============================================================================
// Diagnostic Tile Decoder - Simplified version to trace state machine issues
//==============================================================================

`timescale 1ns / 1ps

module diag_tile_decoder #(
    parameter MAX_WIDTH   = 64,
    parameter MAX_HEIGHT  = 64,
    parameter PIXEL_WIDTH = 10
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     start,
    input  wire [15:0]              frame_width,
    input  wire [15:0]              frame_height,
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
reg [15:0] sb_row, sb_col;
reg [15:0] sb_rows, sb_cols;
reg [15:0] write_offset;
reg [15:0] total_pixels;

// Counters for sub-modules
reg [7:0] entropy_counter;
reg [7:0] itx_counter;
reg [7:0] pred_counter;

integer cycle;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        cycle <= 0;
        sb_row <= 0;
        sb_col <= 0;
        tile_done <= 0;
        write_offset <= 0;
        entropy_counter <= 0;
        itx_counter <= 0;
        pred_counter <= 0;
    end else begin
        cycle <= cycle + 1;
        if (cycle < 20)
            $display("[TIME %0t] cycle=%0d, state=%0d", $time, cycle, state);
        case (state)
            IDLE: begin
                if (cycle < 50)
                    $display("[TIME %0t] IDLE: start=%0b", $time, start);
                tile_done <= 0;
                write_offset <= 0;
                entropy_counter <= 0;
                itx_counter <= 0;
                pred_counter <= 0;
                
                if (start) begin
                    sb_rows <= (frame_height + 63) / 64;
                    sb_cols <= (frame_width + 63) / 64;
                    sb_row <= 0;
                    sb_col <= 0;
                    total_pixels <= frame_width * frame_height;
                    state <= PARSE_SB_HEADER;
                    $display("[TIME %0t] IDLE -> PARSE_SB_HEADER", $time);
                end
            end
            
            PARSE_SB_HEADER: begin
                $display("[TIME %0t] PARSE_SB_HEADER: sb_row=%0d, sb_col=%0d", $time, sb_row, sb_col);
                state <= ENTROPY_DECODE;
            end
            
            ENTROPY_DECODE: begin
                // Simulate entropy decode with counter
                entropy_counter <= entropy_counter + 1;
                if (entropy_counter >= 10) begin
                    entropy_counter <= 0;
                    state <= INVERSE_TX;
                    $display("[TIME %0t] ENTROPY_DECODE -> INVERSE_TX", $time);
                end
            end
            
            INVERSE_TX: begin
                // Simulate inverse transform
                itx_counter <= itx_counter + 1;
                if (itx_counter >= 5) begin
                    itx_counter <= 0;
                    state <= PREDICTION;
                    $display("[TIME %0t] INVERSE_TX -> PREDICTION", $time);
                end
            end
            
            PREDICTION: begin
                // Simulate prediction
                pred_counter <= pred_counter + 1;
                if (pred_counter >= 5) begin
                    pred_counter <= 0;
                    state <= RECONSTRUCTION;
                    $display("[TIME %0t] PREDICTION -> RECONSTRUCTION", $time);
                end
            end
            
            RECONSTRUCTION: begin
                $display("[TIME %0t] RECONSTRUCTION -> CHECK_SB_COMPLETE", $time);
                state <= CHECK_SB_COMPLETE;
            end
            
            CHECK_SB_COMPLETE: begin
                if (sb_col < sb_cols - 1) begin
                    sb_col <= sb_col + 1;
                    state <= PARSE_SB_HEADER;
                    $display("[TIME %0t] CHECK_SB_COMPLETE -> PARSE_SB_HEADER (next col)", $time);
                end else if (sb_row < sb_rows - 1) begin
                    sb_col <= 0;
                    sb_row <= sb_row + 1;
                    state <= PARSE_SB_HEADER;
                    $display("[TIME %0t] CHECK_SB_COMPLETE -> PARSE_SB_HEADER (next row)", $time);
                end else begin
                    state <= WRITE_OUTPUT;
                    $display("[TIME %0t] CHECK_SB_COMPLETE -> WRITE_OUTPUT", $time);
                end
            end
            
            WRITE_OUTPUT: begin
                if (write_offset < total_pixels) begin
                    write_offset <= write_offset + 16;
                    if (write_offset % 256 == 0)
                        $display("[TIME %0t] WRITE_OUTPUT: offset=%0d / %0d", $time, write_offset, total_pixels);
                end else begin
                    state <= DONE;
                    $display("[TIME %0t] WRITE_OUTPUT -> DONE", $time);
                end
            end
            
            DONE: begin
                tile_done <= 1;
                $display("[TIME %0t] DONE! Total cycles=%0d", $time, cycle);
                state <= IDLE;
            end
        endcase
    end
end

endmodule
