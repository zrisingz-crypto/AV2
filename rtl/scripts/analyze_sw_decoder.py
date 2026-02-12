#!/usr/bin/env python3
"""
Analyze software decoder output to understand YUV format and data
"""

import struct
import sys

def read_yuv420(filename, width, height):
    """
    Read YUV420 file
    
    YUV420 format:
    - Y plane: width * height bytes
    - U plane: (width/2) * (height/2) bytes
    - V plane: (width/2) * (height/2) bytes
    """
    with open(filename, 'rb') as f:
        data = f.read()
    
    y_size = width * height
    uv_size = (width // 2) * (height // 2)
    
    if len(data) < y_size + 2 * uv_size:
        print(f"Error: File too small. Expected at least {y_size + 2 * uv_size} bytes, got {len(data)}")
        sys.exit(1)
    
    y_plane = data[:y_size]
    u_plane = data[y_size:y_size + uv_size]
    v_plane = data[y_size + uv_size:y_size + 2 * uv_size]
    
    return y_plane, u_plane, v_plane

def print_stats(plane, name):
    """Print statistics for a plane"""
    if len(plane) == 0:
        return
    
    pixels = [b for b in plane]
    
    print(f"\n{name} Plane Statistics:")
    print(f"  Size: {len(pixels)} bytes")
    print(f"  Min: {min(pixels)}")
    print(f"  Max: {max(pixels)}")
    print(f"  Mean: {sum(pixels) / len(pixels):.2f}")
    print(f"  First 10 pixels: {pixels[:10]}")
    print(f"  Last 10 pixels: {pixels[-10:]}")
    
    # Count unique values
    unique = set(pixels)
    print(f"  Unique values: {len(unique)}")
    if len(unique) <= 256:
        print(f"  Value range: {min(unique)} - {max(unique)}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python analyze_sw_decoder.py <yuv_file> <width> <height>")
        sys.exit(1)
    
    yuv_file = sys.argv[1]
    width = int(sys.argv[2])
    height = int(sys.argv[3])
    
    print(f"Analyzing: {yuv_file}")
    print(f"Resolution: {width}x{height}")
    print(f"Expected YUV420 size: {width * height * 3 // 2} bytes")
    
    # Read YUV420 file
    y_plane, u_plane, v_plane = read_yuv420(yuv_file, width, height)
    
    # Print statistics
    print_stats(y_plane, "Y")
    print_stats(u_plane, "U")
    print_stats(v_plane, "V")
    
    # Calculate expected sizes
    y_size = width * height
    uv_size = (width // 2) * (height // 2)
    
    print(f"\nExpected YUV420 format:")
    print(f"  Y plane: {y_size} bytes ({width}x{height})")
    print(f"  U plane: {uv_size} bytes ({width//2}x{height//2})")
    print(f"  V plane: {uv_size} bytes ({width//2}x{height//2})")
    print(f"  Total: {y_size + 2 * uv_size} bytes")
    
    # Check if this is a simple gradient or has actual image data
    check_len = min(256, len(y_plane))
    is_gradient = all(y_plane[i] == i % 256 for i in range(check_len))
    
    if is_gradient:
        print("\nNote: Y plane appears to be a simple gradient (0,1,2,3,...)")
        print("This is test data, not real image content.")
    else:
        print("\nNote: Y plane contains actual image data (not a simple gradient)")

if __name__ == "__main__":
    main()