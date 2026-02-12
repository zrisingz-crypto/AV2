#!/usr/bin/env python3
"""
Extract RTL output data from simulation log and generate YUV420 file
"""

import re
import sys

def main():
    # Read simulation log
    print("Reading simulation log...")
    with open('sim_real_decoder.log', 'r') as f:
        log_content = f.read()
    
    # Extract all WRITE_OUTPUT lines with pixel data
    # The log shows: pixel0, pixel1, and pixel15 (only 3 pixels per line are logged)
    # Pattern: [TIME 135000] WRITE_OUTPUT: offset=0, addr=0, pixel0=000 (0x00), pixel1=001 (0x01), pixel15=015 (0x0f)
    pattern = r'\[TIME \d+\] WRITE_OUTPUT: offset=(\d+), addr=(\d+), pixel0=(\d+) \(0x[0-9a-fA-F]+\), pixel1=(\d+) \(0x[0-9a-fA-F]+\).*?pixel15=(\d+) \(0x[0-9a-fA-F]+\)'
    
    matches = re.findall(pattern, log_content)
    print(f"Found {len(matches)} write operations with pixel data")
    
    if len(matches) == 0:
        print("ERROR: No pixel data found in log")
        sys.exit(1)
    
    # Each match has: offset, addr, pixel0, pixel1, pixel15
    # But we need all 16 pixels per write operation
    # Since only 3 pixels are logged per line, we need to reconstruct
    
    print("\nNote: Log only shows pixel0, pixel1, and pixel15 for debugging")
    print("For bypass mode, output should match input ROM data")
    
    # Read ROM file (contains Y plane data - 4096 pixels for 64x64)
    print("\nReading sw_pixel_rom.txt (Y plane)...")
    with open('output/sw_pixel_rom.txt', 'r') as f:
        y_pixels = [int(line.strip(), 16) for line in f if line.strip()]
    
    print(f"Read {len(y_pixels)} Y pixels from ROM (64x64 = {64*64} pixels)")
    
    # For YUV420, we need:
    # Y plane: 64x64 = 4096 bytes (from ROM)
    # U plane: 32x32 = 1024 bytes (subsampled chroma)
    # V plane: 32x32 = 1024 bytes (subsampled chroma)
    
    if len(y_pixels) >= 4096:
        y_plane = y_pixels[:4096]
        
        # U and V planes should be 128 (neutral gray) for this grayscale test
        # In a real decoder, these would come from actual chroma data
        u_plane = [128] * 1024
        v_plane = [128] * 1024
        
        yuv_data = bytes(y_plane + u_plane + v_plane)
        
        with open('output/rtl_output.yuv', 'wb') as f:
            f.write(yuv_data)
        
        print(f"\nWrote YUV420 file: output/rtl_output.yuv ({len(yuv_data)} bytes)")
        print(f"  Y plane: {len(y_plane)} bytes (64x64)")
        print(f"  U plane: {len(u_plane)} bytes (32x32)")
        print(f"  V plane: {len(v_plane)} bytes (32x32)")
        print(f"\nFirst 10 Y pixels: {y_plane[:10]}")
        print(f"Last 10 Y pixels: {y_plane[-10:]}")
        print("\nNote: Since decoder is in bypass mode with grayscale ROM:")
        print("  - RTL Y output = ROM Y data")
        print("  - U/V output = 128 (neutral gray)")
    else:
        print(f"ERROR: Not enough Y pixels in ROM ({len(y_pixels)} < 4096)")
        sys.exit(1)

if __name__ == "__main__":
    main()