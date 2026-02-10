#!/usr/bin/env python3
"""Convert RTL text output to YUV binary format"""

import sys
import struct

def convert_rtl_to_yuv(rtl_file, yuv_file, width=64, height=64):
    # Initialize Y plane with zeros
    y_plane = [0] * (width * height)
    
    # Read RTL output
    with open(rtl_file, 'r') as f:
        for line in f:
            line = line.strip()
            if 'pixel[' in line and '] = ' in line:
                # Parse: pixel[addr][offset] = value
                try:
                    parts = line.split(' = ')
                    if len(parts) == 2:
                        value_str = parts[1]
                        if value_str == 'x':
                            value = 0  # Skip undefined values
                        else:
                            value = int(value_str)
                        
                        # Parse address
                        addr_part = parts[0]
                        # Extract addr and offset from pixel[addr][offset]
                        import re
                        match = re.match(r'pixel\[(\d+)\]\[(\d+)\]', addr_part)
                        if match:
                            addr = int(match.group(1))
                            offset = int(match.group(2))
                            # Calculate pixel index
                            pixel_idx = addr + offset
                            if 0 <= pixel_idx < width * height:
                                y_plane[pixel_idx] = value & 0xFF
                except Exception as e:
                    print(f"Error parsing line: {line}, error: {e}")
    
    # Write YUV file (Y only for now, U and V are 128)
    with open(yuv_file, 'wb') as f:
        # Y plane
        for val in y_plane:
            f.write(struct.pack('B', val))
        # U plane (all 128)
        for _ in range(width * height // 4):
            f.write(struct.pack('B', 128))
        # V plane (all 128)
        for _ in range(width * height // 4):
            f.write(struct.pack('B', 128))
    
    print(f"Converted {rtl_file} -> {yuv_file}")
    print(f"Y plane: {len(y_plane)} pixels")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: convert_rtl_to_yuv.py <rtl_output.txt> <output.yuv>")
        sys.exit(1)
    
    convert_rtl_to_yuv(sys.argv[1], sys.argv[2])
