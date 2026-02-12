#!/usr/bin/env python3
"""Compare Y plane of two YUV420 files for AV2 test."""

import argparse

def compare_yuv_files(file1, file2, width=64, height=64):
    """Compare Y plane of two YUV files."""
    y_size = width * height
    
    with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
        y1 = f1.read(y_size)
        y2 = f2.read(y_size)
        
        if len(y1) < y_size or len(y2) < y_size:
            print(f"Error: Files too small. Need at least {y_size} bytes for Y plane")
            return False
        
        mismatches = 0
        for i in range(y_size):
            if y1[i] != y2[i]:
                mismatches += 1
        
        total = y_size
        match_percent = (1 - mismatches / total) * 100
        
        print(f"============================================")
        print(f"Y Plane Comparison Results ({width}x{height})")
        print(f"============================================")
        print(f"File 1: {file1}")
        print(f"File 2: {file2}")
        print(f"Total pixels: {total}")
        print(f"Mismatched pixels: {mismatches} ({100 - match_percent:.2f}%)")
        print(f"Match rate: {match_percent:.2f}%")
        print(f"============================================")
        
        return mismatches == 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Y plane of two YUV420 files")
    parser.add_argument("file1", help="First YUV file (RTL output)")
    parser.add_argument("file2", help="Second YUV file (SW output)")
    parser.add_argument("-w", "--width", type=int, default=64, help="Frame width (default: 64)")
    parser.add_argument("-H", "--height", type=int, default=64, help="Frame height (default: 64)")
    
    args = parser.parse_args()
    
    if compare_yuv_files(args.file1, args.file2, args.width, args.height):
        print("\n✅ Y planes match perfectly!")
    else:
        print("\n❌ Y planes do NOT match!")
