#!/usr/bin/env python3
"""
Compare RTL decoder output with software decoder output
"""

import struct

def read_yuv420(filename, width=64, height=64):
    """Read YUV420 format file"""
    with open(filename, 'rb') as f:
        data = f.read()
    
    y_size = width * height
    uv_size = y_size // 4
    
    y_plane = data[:y_size]
    u_plane = data[y_size:y_size + uv_size]
    v_plane = data[y_size + uv_size:y_size + 2*uv_size]
    
    return y_plane, u_plane, v_plane

def compare_plane(soft_plane, hw_plane, plane_name):
    """Compare two planes and return statistics"""
    if len(soft_plane) != len(hw_plane):
        print(f"  {plane_name} plane size mismatch: soft={len(soft_plane)}, hw={len(hw_plane)}")
        return None
    
    total_pixels = len(soft_plane)
    different_pixels = 0
    total_diff = 0
    max_diff = 0
    
    for i in range(total_pixels):
        diff = abs(soft_plane[i] - hw_plane[i])
        if diff > 0:
            different_pixels += 1
            total_diff += diff
            if diff > max_diff:
                max_diff = diff
    
    return {
        'total_pixels': total_pixels,
        'different_pixels': different_pixels,
        'percentage': (different_pixels / total_pixels * 100) if total_pixels > 0 else 0,
        'avg_diff': total_diff / different_pixels if different_pixels > 0 else 0,
        'max_diff': max_diff
    }

def main():
    print("=" * 60)
    print("RTL Decoder vs Software Decoder Comparison")
    print("=" * 60)
    
    # Read software decoder output
    print("\n[Software Decoder]")
    try:
        soft_y, soft_u, soft_v = read_yuv420('sw_decoder_out_0.yuv')
        print(f"  ✓ Read sw_decoder_out_0.yuv")
        print(f"    Y plane: {len(soft_y)} bytes")
        print(f"    U plane: {len(soft_u)} bytes")
        print(f"    V plane: {len(soft_v)} bytes")
    except Exception as e:
        print(f"  ✗ Failed to read software decoder output: {e}")
        return
    
    # Read RTL decoder output
    print("\n[RTL Decoder]")
    try:
        hw_y, hw_u, hw_v = read_yuv420('hw_decoder_out_0_fixed.yuv')
        print(f"  ✓ Read hw_decoder_out_0_fixed.yuv")
        print(f"    Y plane: {len(hw_y)} bytes")
        print(f"    U plane: {len(hw_u)} bytes")
        print(f"    V plane: {len(hw_v)} bytes")
    except Exception as e:
        print(f"  ✗ Failed to read RTL decoder output: {e}")
        return
    
    # Compare Y plane
    print("\n[Y Plane Comparison]")
    y_stats = compare_plane(soft_y, hw_y, 'Y')
    if y_stats:
        print(f"  Total pixels: {y_stats['total_pixels']}")
        print(f"  Different pixels: {y_stats['different_pixels']} ({y_stats['percentage']:.2f}%)")
        print(f"  Average difference: {y_stats['avg_diff']:.2f}")
        print(f"  Maximum difference: {y_stats['max_diff']}")
    
    # Compare U plane
    print("\n[U Plane Comparison]")
    u_stats = compare_plane(soft_u, hw_u, 'U')
    if u_stats:
        print(f"  Total pixels: {u_stats['total_pixels']}")
        print(f"  Different pixels: {u_stats['different_pixels']} ({u_stats['percentage']:.2f}%)")
        print(f"  Average difference: {u_stats['avg_diff']:.2f}")
        print(f"  Maximum difference: {u_stats['max_diff']}")
    
    # Compare V plane
    print("\n[V Plane Comparison]")
    v_stats = compare_plane(soft_v, hw_v, 'V')
    if v_stats:
        print(f"  Total pixels: {v_stats['total_pixels']}")
        print(f"  Different pixels: {v_stats['different_pixels']} ({v_stats['percentage']:.2f}%)")
        print(f"  Average difference: {v_stats['avg_diff']:.2f}")
        print(f"  Maximum difference: {v_stats['max_diff']}")
    
    # Overall summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    if y_stats and u_stats and v_stats:
        total_diff_pixels = y_stats['different_pixels'] + u_stats['different_pixels'] + v_stats['different_pixels']
        total_pixels = y_stats['total_pixels'] + u_stats['total_pixels'] + v_stats['total_pixels']
        match_percentage = ((total_pixels - total_diff_pixels) / total_pixels * 100) if total_pixels > 0 else 0
        
        print(f"  Total pixels: {total_pixels}")
        print(f"  Matching pixels: {total_pixels - total_diff_pixels} ({match_percentage:.2f}%)")
        print(f"  Different pixels: {total_diff_pixels} ({100 - match_percentage:.2f}%)")
        
        if match_percentage > 95:
            print(f"\n  ✗ RESULT: POOR MATCH ({match_percentage:.2f}%)")
            print(f"     Note: RTL uses stub tile decoder (outputs zeros)")
        elif match_percentage > 80:
            print(f"\n  ✗ RESULT: POOR MATCH ({match_percentage:.2f}%)")
        else:
            print(f"\n  ✗ RESULT: POOR MATCH ({match_percentage:.2f}%)")
    
    print("=" * 60)

if __name__ == '__main__':
    main()