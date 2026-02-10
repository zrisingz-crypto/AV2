#!/usr/bin/env python3
"""Compare two YUV files"""

import sys
import struct

def read_yuv(yuv_file, width=64, height=64):
    with open(yuv_file, 'rb') as f:
        data = f.read()
    
    y_size = width * height
    uv_size = y_size // 4
    
    y_plane = list(data[0:y_size])
    u_plane = list(data[y_size:y_size+uv_size])
    v_plane = list(data[y_size+uv_size:y_size+2*uv_size])
    
    return y_plane, u_plane, v_plane

def compare_yuv(file1, file2, width=64, height=64):
    y1, u1, v1 = read_yuv(file1, width, height)
    y2, u2, v2 = read_yuv(file2, width, height)
    
    print(f"Comparing {file1} vs {file2}")
    print(f"Resolution: {width}x{height}")
    print()
    
    # Compare Y plane
    diff_count = 0
    max_diff = 0
    first_diff = None
    
    for i in range(len(y1)):
        diff = abs(y1[i] - y2[i])
        if diff > 0:
            diff_count += 1
            if diff > max_diff:
                max_diff = diff
            if first_diff is None:
                first_diff = (i, y1[i], y2[i])
    
    print(f"Y plane comparison:")
    print(f"  Total pixels: {len(y1)}")
    print(f"  Different pixels: {diff_count} ({100*diff_count/len(y1):.2f}%)")
    print(f"  Max difference: {max_diff}")
    if first_diff:
        x = first_diff[0] % width
        y = first_diff[0] // width
        print(f"  First diff at ({x}, {y}): {file1}={first_diff[1]}, {file2}={first_diff[2]}")
    print()
    
    # Print first 64 Y values from each file
    print(f"First 64 Y values:")
    print(f"  {file1}: {y1[:64]}")
    print(f"  {file2}: {y2[:64]}")
    
    return diff_count == 0

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: compare_yuv.py <file1.yuv> <file2.yuv>")
        sys.exit(1)
    
    match = compare_yuv(sys.argv[1], sys.argv[2])
    print()
    if match:
        print("✓ Files match!")
    else:
        print("✗ Files differ")
    sys.exit(0 if match else 1)
