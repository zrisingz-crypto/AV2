#!/usr/bin/env python3
"""
Compare RTL and SW decoder outputs
"""

import sys
import os

def read_yuv_file(filename, width, height, num_frames=2):
    """Read YUV file and return list of frames"""
    y_size = width * height
    uv_size = y_size // 4
    frame_size = y_size + 2 * uv_size  # YUV420
    
    with open(filename, 'rb') as f:
        data = f.read()
    
    frames = []
    for i in range(num_frames):
        start = i * frame_size
        end = start + frame_size
        if end <= len(data):
            frame = data[start:end]
            frames.append(frame)
    
    return frames

def compare_frames(rtl_frame, sw_frame, width, height):
    """Compare two frames and return statistics"""
    y_size = width * height
    
    # Compare Y plane
    y_diffs = []
    for i in range(y_size):
        diff = abs(rtl_frame[i] - sw_frame[i])
        if diff > 0:
            y_diffs.append((i, diff))
    
    # Compare U plane
    u_start = y_size
    u_diffs = []
    for i in range(y_size // 4):
        diff = abs(rtl_frame[u_start + i] - sw_frame[u_start + i])
        if diff > 0:
            u_diffs.append((i, diff))
    
    # Compare V plane
    v_start = y_size + (y_size // 4)
    v_diffs = []
    for i in range(y_size // 4):
        diff = abs(rtl_frame[v_start + i] - sw_frame[v_start + i])
        if diff > 0:
            v_diffs.append((i, diff))
    
    return {
        'y_diffs': y_diffs,
        'u_diffs': u_diffs,
        'v_diffs': v_diffs,
        'total_diffs': len(y_diffs) + len(u_diffs) + len(v_diffs)
    }

def main():
    width = 64
    height = 64
    
    # Read RTL output
    print("Reading RTL decoder output...")
    rtl_frames = read_yuv_file('hw_decoder_out_0.yuv', width, height, num_frames=1)
    
    # Read SW output
    print("Reading SW decoder output...")
    sw_frames = read_yuv_file('sw_decoder_out_0.yuv', width, height, num_frames=1)
    
    if len(rtl_frames) == 0:
        print("Error: No frames found in RTL output")
        return 1
    
    if len(sw_frames) == 0:
        print("Error: No frames found in SW output")
        return 1
    
    print(f"\nRTL frames: {len(rtl_frames)}")
    print(f"SW frames: {len(sw_frames)}")
    
    # Compare first frame
    print("\n" + "="*60)
    print(f"Comparing Frame 0 ({width}x{height})")
    print("="*60)
    
    result = compare_frames(rtl_frames[0], sw_frames[0], width, height)
    
    y_total = width * height
    u_total = y_total // 4
    v_total = y_total // 4
    
    print(f"\nY plane: {len(result['y_diffs'])}/{y_total} pixels differ")
    print(f"U plane: {len(result['u_diffs'])}/{u_total} pixels differ")
    print(f"V plane: {len(result['v_diffs'])}/{v_total} pixels differ")
    print(f"Total: {result['total_diffs']}/{y_total + u_total + v_total} pixels differ")
    
    if result['total_diffs'] > 0:
        print(f"\nSample differences (first 20):")
        
        # Show some Y differences
        if result['y_diffs']:
            print("\n  Y plane differences:")
            for idx, diff in result['y_diffs'][:10]:
                row = idx // width
                col = idx % width
                rtl_val = rtl_frames[0][idx]
                sw_val = sw_frames[0][idx]
                print(f"    [{row:2d},{col:2d}] RTL={rtl_val:3d} SW={sw_val:3d} diff={diff:3d}")
        
        # Show some U differences
        if result['u_diffs']:
            print("\n  U plane differences:")
            for idx, diff in result['u_diffs'][:10]:
                row = idx // (width // 2)
                col = idx % (width // 2)
                rtl_val = rtl_frames[0][y_total + idx]
                sw_val = sw_frames[0][y_total + idx]
                print(f"    [{row:2d},{col:2d}] RTL={rtl_val:3d} SW={sw_val:3d} diff={diff:3d}")
        
        # Show some V differences
        if result['v_diffs']:
            print("\n  V plane differences:")
            for idx, diff in result['v_diffs'][:10]:
                row = idx // (width // 2)
                col = idx % (width // 2)
                rtl_val = rtl_frames[0][y_total + (y_total//4) + idx]
                sw_val = sw_frames[0][y_total + (y_total//4) + idx]
                print(f"    [{row:2d},{col:2d}] RTL={rtl_val:3d} SW={sw_val:3d} diff={diff:3d}")
    
    # Calculate PSNR
    if result['total_diffs'] > 0:
        mse = 0
        for i in range(len(rtl_frames[0])):
            diff = rtl_frames[0][i] - sw_frames[0][i]
            mse += diff * diff
        mse /= len(rtl_frames[0])
        psnr = 10 * (10 * (255*255) / mse) if mse > 0 else float('inf')
        print(f"\nMSE: {mse:.2f}")
        print(f"PSNR: {psnr:.2f} dB")
    
    print("\n" + "="*60)
    
    if result['total_diffs'] == 0:
        print("✓ Frames match perfectly!")
        return 0
    else:
        print(f"✗ Frames differ by {result['total_diffs']} pixels")
        return 1

if __name__ == '__main__':
    sys.exit(main())