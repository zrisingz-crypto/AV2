#!/usr/bin/env python3
"""
Simple YUV file comparison script
Compares two YUV420 files pixel by pixel
"""

import sys

def read_yuv420(filename, width, height):
    """Read YUV420 file and return Y, U, V planes"""
    y_size = width * height
    uv_size = width * height // 4
    
    with open(filename, 'rb') as f:
        data = f.read()
    
    # Y plane (full resolution)
    y = list(data[:y_size])
    
    # U plane (1/4 resolution)
    u = list(data[y_size:y_size + uv_size])
    
    # V plane (1/4 resolution)
    v = list(data[y_size + uv_size:y_size + 2*uv_size])
    
    return y, u, v

def compare_pixels(ref, test, name):
    """Compare two pixel arrays and return mismatch percentage"""
    if len(ref) != len(test):
        print(f"ERROR: {name} plane size mismatch: ref={len(ref)}, test={len(test)}")
        return 100.0
    
    mismatches = 0
    for i in range(len(ref)):
        diff = abs(ref[i] - test[i])
        if diff > 5:  # Allow small tolerance
            mismatches += 1
    
    mismatch_pct = (mismatches / len(ref)) * 100
    return mismatch_pct

def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <ref_yuv> <test_yuv> <width> <height>")
        sys.exit(1)
    
    ref_file = sys.argv[1]
    test_file = sys.argv[2]
    width = int(sys.argv[3])
    height = int(sys.argv[4])
    
    print(f"Comparing {ref_file} vs {test_file}")
    print(f"Resolution: {width}x{height}")
    print()
    
    # Read files
    ref_y, ref_u, ref_v = read_yuv420(ref_file, width, height)
    test_y, test_u, test_v = read_yuv420(test_file, width, height)
    
    # Compare each plane
    y_mismatch = compare_pixels(ref_y, test_y, "Y")
    u_mismatch = compare_pixels(ref_u, test_u, "U")
    v_mismatch = compare_pixels(ref_v, test_v, "V")
    
    # Calculate overall average
    total_mismatch = (y_mismatch + u_mismatch + v_mismatch) / 3.0
    
    # Print results
    print("Plane Comparison Results:")
    print(f"  Y plane: {y_mismatch:.2f}% mismatch")
    print(f"  U plane: {u_mismatch:.2f}% mismatch")
    print(f"  V plane: {v_mismatch:.2f}% mismatch")
    print(f"  Overall: {total_mismatch:.2f}% mismatch")
    print()
    
    # Verdict
    if total_mismatch < 1.0:
        print("✓ EXCELLENT: Output matches reference very well (< 1% mismatch)")
    elif total_mismatch < 5.0:
        print("✓ GOOD: Output matches reference well (< 5% mismatch)")
    elif total_mismatch < 10.0:
        print("⚠ ACCEPTABLE: Output has moderate mismatch (< 10% mismatch)")
    else:
        print("✗ POOR: Output has significant mismatch (> 10% mismatch)")
    
    return total_mismatch

if __name__ == "__main__":
    main()