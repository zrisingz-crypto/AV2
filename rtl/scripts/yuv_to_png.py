#!/usr/bin/env python3
"""
Convert YUV420 raw file to PNG image (without numpy dependency)
Supports 64x64 resolution (first frame only)
"""

import sys
import struct
from PIL import Image

def yuv420_to_png(yuv_file, png_file, width=64, height=64):
    """Convert YUV420 raw file to PNG"""
    
    # Calculate plane sizes
    y_size = width * height
    uv_size = (width // 2) * (height // 2)
    
    # Read YUV file
    with open(yuv_file, 'rb') as f:
        data = f.read()
    
    if len(data) < y_size + 2 * uv_size:
        print(f"Warning: File smaller than expected YUV420 format")
        print(f"Expected: {y_size + 2*uv_size} bytes, Got: {len(data)} bytes")
        if len(data) < y_size:
            print("Error: Not enough data for Y plane")
            return None
    
    # Extract Y plane
    y_data = list(data[:y_size])
    
    # Extract U and V planes
    u_data = list(data[y_size:y_size + uv_size])
    v_data = list(data[y_size + uv_size:y_size + 2 * uv_size])
    
    # If U/V planes are missing or incomplete, fill with 128
    if len(u_data) < uv_size:
        u_data = [128] * uv_size
    if len(v_data) < uv_size:
        v_data = [128] * uv_size
    
    # Create RGB image
    rgb_data = []
    for y in range(height):
        for x in range(width):
            # Get Y value
            Y = y_data[y * width + x]
            
            # Get U and V (downsampled by 2x2)
            u_idx = (y // 2) * (width // 2) + (x // 2)
            U = u_data[u_idx] - 128
            V = v_data[u_idx] - 128
            
            # YUV to RGB conversion (BT.601)
            R = int(Y + 1.402 * V)
            G = int(Y - 0.344136 * U - 0.714136 * V)
            B = int(Y + 1.772 * U)
            
            # Clip to 0-255
            R = max(0, min(255, R))
            G = max(0, min(255, G))
            B = max(0, min(255, B))
            
            rgb_data.append((R, G, B))
    
    # Save as PNG
    img = Image.new('RGB', (width, height))
    img.putdata(rgb_data)
    img.save(png_file)
    
    print(f"Converted: {yuv_file} -> {png_file}")
    print(f"Resolution: {width}x{height}")
    
    return img

def y_only_to_png(yuv_file, png_file, width=64, height=64):
    """Convert Y-only (grayscale) raw file to PNG"""
    
    # Read Y plane only
    with open(yuv_file, 'rb') as f:
        y_data = f.read(width * height)
    
    if len(y_data) < width * height:
        print(f"Error: File too small. Expected {width*height} bytes, got {len(y_data)}")
        return None
    
    # Create grayscale image
    img = Image.new('L', (width, height))
    img.putdata(list(y_data[:width*height]))
    img.save(png_file)
    
    print(f"Converted (grayscale): {yuv_file} -> {png_file}")
    print(f"Resolution: {width}x{height}")
    
    return img

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: yuv_to_png.py <input.yuv> <output.png> [width] [height] [--y-only]")
        print("  --y-only: Treat input as Y-only (grayscale) data")
        print("\nExample:")
        print("  python3 yuv_to_png.py hw_decoder_out_real.yuv output.png")
        sys.exit(1)
    
    yuv_file = sys.argv[1]
    png_file = sys.argv[2]
    width = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    height = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    y_only = "--y-only" in sys.argv
    
    if y_only:
        y_only_to_png(yuv_file, png_file, width, height)
    else:
        yuv420_to_png(yuv_file, png_file, width, height)
