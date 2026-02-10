//==============================================================================
// Output Controller Module - Simple and Reliable
//==============================================================================

`timescale 1ns / 1ps

module av2_output_ctrl #(
    parameter DATA_WIDTH = 128
)(
    input  wire                    clk,
    input  wire                    rst_n,
    
    input  wire                    start,
    input  wire [15:0]             frame_width,
    input  wire [15:0]             frame_height,
    output wire [31:0]             fb_rd_addr,
    output wire [1:0]              fb_rd_sel_plane,  // 0: Y, 1: U, 2: V
    input  wire [DATA_WIDTH-1:0]   fb_rd_data,
    
    // AXI4-Stream Output
    output wire [DATA_WIDTH-1:0]   m_axis_tdata,
    output reg                     m_axis_tvalid,
    input  wire                    m_axis_tready,
    output reg                     m_axis_tlast,
    output wire [15:0]             m_axis_tkeep,
    output reg  [1:0]              m_axis_tuser,
    
    // Status output
    output wire                    done
);

// Calculate pixel counts for YUV420 planes
wire [15:0] y_pixels = frame_width * frame_height;           // 4096 for 64x64
wire [15:0] uv_pixels = (frame_width * frame_height) >> 2;   // 1024 for 64x64
wire [15:0] total_pixels = y_pixels + (uv_pixels << 1);      // 6144 for 64x64

// Determine current plane based on global counter
wire [1:0] current_plane;
wire [15:0] pixels_in_current_plane;
wire [15:0] offset_in_current_plane;

assign current_plane = (global_count < y_pixels) ? 2'd0 :
                      (global_count < y_pixels + uv_pixels) ? 2'd1 : 2'd2;

assign pixels_in_current_plane = (current_plane == 2'd0) ? y_pixels : uv_pixels;

assign offset_in_current_plane = (current_plane == 2'd0) ? global_count :
                                 (current_plane == 2'd1) ? (global_count - y_pixels) :
                                 (global_count - y_pixels - uv_pixels);

// Global counter for all pixels
reg [15:0] global_count;
reg [31:0] rd_addr;
reg transfer_active;
reg final_cycle; // Flag to indicate we sent tlast

// Calculate remaining pixels in current plane
wire [15:0] remaining_in_plane = pixels_in_current_plane - offset_in_current_plane;

// Calculate how many pixels to transfer (max 16)
wire [4:0] pixels_to_transfer = (remaining_in_plane >= 16) ? 5'd16 : remaining_in_plane[4:0];

// tkeep based on pixels to transfer
wire [15:0] tkeep_wire;
assign tkeep_wire = (pixels_to_transfer == 5'd1)  ? 16'h0001 :
                   (pixels_to_transfer == 5'd2)  ? 16'h0003 :
                   (pixels_to_transfer == 5'd3)  ? 16'h0007 :
                   (pixels_to_transfer == 5'd4)  ? 16'h000F :
                   (pixels_to_transfer == 5'd5)  ? 16'h001F :
                   (pixels_to_transfer == 5'd6)  ? 16'h003F :
                   (pixels_to_transfer == 5'd7)  ? 16'h007F :
                   (pixels_to_transfer == 5'd8)  ? 16'h00FF :
                   (pixels_to_transfer == 5'd9)  ? 16'h01FF :
                   (pixels_to_transfer == 5'd10) ? 16'h03FF :
                   (pixels_to_transfer == 5'd11) ? 16'h07FF :
                   (pixels_to_transfer == 5'd12) ? 16'h0FFF :
                   (pixels_to_transfer == 5'd13) ? 16'h1FFF :
                   (pixels_to_transfer == 5'd14) ? 16'h3FFF :
                   (pixels_to_transfer == 5'd15) ? 16'h7FFF :
                   (pixels_to_transfer == 5'd16) ? 16'hFFFF : 16'h0000;

// Write cycle counter for debugging
reg [15:0] write_cycle_count;

reg done_r;
reg start_r;  // Previous value of start for edge detection

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        transfer_active <= 1'b0;
        global_count <= 16'd0;
        rd_addr <= 32'd0;
        m_axis_tvalid <= 1'b0;
        m_axis_tlast <= 1'b0;
        m_axis_tuser <= 2'b00;
        final_cycle <= 1'b0;
        done_r <= 1'b0;
        start_r <= 1'b0;
        write_cycle_count <= 16'd0;
    end else begin
        done_r <= 1'b0; // Default: done is only a pulse
        start_r <= start;  // Register start signal
        
        // Debug: detect start edge
        if (start && !start_r) begin
            $display("[%0t] Output Ctrl: START EDGE detected! start=%d, start_r=%d, transfer_active=%d, done_r=%d", 
                     $time, start, start_r, transfer_active, done_r);
        end
        
        // Debug: print every cycle during active transfer
        if (transfer_active) begin
            $display("[%0t] Output Ctrl: active=%d, valid=%d, ready=%d, count=%0d/%0d, addr=%0d, plane=%0d, remaining=%0d, pixels_xfer=%0d", 
                     $time, transfer_active, m_axis_tvalid, m_axis_tready, global_count, total_pixels, rd_addr, current_plane, remaining_in_plane, pixels_to_transfer);
            if (m_axis_tvalid && m_axis_tready) begin
                $display("[%0t] Output Ctrl:   transfer_handshake - will_count=%0d, final=%0d, offset=%0d", 
                         $time, global_count + pixels_to_transfer, final_cycle, offset_in_current_plane);
            end
        end
        
        // Debug: periodically print state even when not active
        if (!transfer_active && ($time % 1024 == 0)) begin
            $display("[%0t] Output Ctrl: IDLE - start=%d, transfer_active=%d, done_r=%d", 
                     $time, start, transfer_active, done_r);
        end
        
        // Start signal initiates frame transfer (only if not already active)
        if (start && !transfer_active && !done_r) begin
            transfer_active <= 1'b1;
            global_count <= 16'd0;
            rd_addr <= 32'd0;
            m_axis_tvalid <= 1'b1;
            m_axis_tuser <= 2'b10; // frame_start
            final_cycle <= 1'b0;
            write_cycle_count <= 16'd0;
            $display("[%0t] Output Ctrl: START detected, total_pixels=%0d", $time, total_pixels);
        end
        
        // Active transfer
        if (transfer_active && m_axis_tvalid && m_axis_tready) begin
            // Increment write cycle counter
            write_cycle_count <= write_cycle_count + 16'd1;
            
            // If we're in the final cycle (after sending tlast), stop
            if (final_cycle) begin
                transfer_active <= 1'b0;
                m_axis_tvalid <= 1'b0;
                m_axis_tlast <= 1'b0;
                m_axis_tuser <= 2'b00;
                final_cycle <= 1'b0;
                done_r <= 1'b1; // Assert done for one cycle
                $display("[%0t] Output Ctrl: Frame complete (with valid) at count=%0d, write_cycles=%0d", 
                         $time, global_count, write_cycle_count);
            end else begin
                // Check if this is the final transfer
                if (global_count + pixels_to_transfer >= total_pixels) begin
                    // Final transfer - set tlast, frame_end, but keep valid=1 for handshake
                    m_axis_tlast <= 1'b1;
                    m_axis_tuser <= 2'b01; // frame_end
                    final_cycle <= 1'b1;
                    // Keep valid=1 and transfer_active=1 to allow handshake
                    // Don't set valid=0 or transfer_active=0 yet
                    done_r <= 1'b0; // Don't set done yet
                    $display("[%0t] Output Ctrl: Final transfer, count=%0d->%0d, total=%0d, write_cycles=%0d",
                             $time, global_count, global_count + pixels_to_transfer, total_pixels, write_cycle_count + 1);
                end else begin
                    // Normal transfer - update counters
                    global_count <= global_count + pixels_to_transfer;
                    
                    // Update address - increment only within each plane
                    if (offset_in_current_plane + pixels_to_transfer >= pixels_in_current_plane) begin
                        // End of current plane, reset address for next plane
                        rd_addr <= 32'd0;
                    end else begin
                        // Still in current plane, increment address
                        rd_addr <= rd_addr + 32'd1;
                    end
                    
                    m_axis_tlast <= 1'b0;
                    m_axis_tuser <= 2'b00;
                end
            end
        end else if (transfer_active) begin
            // Transfer is active but not ready - this includes the final_cycle case
            if (final_cycle) begin
                $display("[%0t] Output Ctrl: Stopping (no valid ready), final_cycle=1", $time);
                transfer_active <= 1'b0;
                m_axis_tvalid <= 1'b0;
                m_axis_tlast <= 1'b0;
                m_axis_tuser <= 2'b00;
                final_cycle <= 1'b0;
                done_r <= 1'b1;
            end else if (!m_axis_tvalid) begin
                $display("[%0t] Output Ctrl: Waiting (not valid)", $time);
            end else begin
                $display("[%0t] Output Ctrl: Waiting (not ready), valid=%d", $time, m_axis_tvalid);
            end
        end else begin
            done_r <= 1'b0;
        end
    end
end

assign fb_rd_addr = rd_addr;

// Map current plane to fb_rd_sel_plane
assign fb_rd_sel_plane = current_plane;

// Output data from current plane
assign m_axis_tdata = fb_rd_data;
assign m_axis_tkeep = tkeep_wire;
assign done = done_r;

endmodule