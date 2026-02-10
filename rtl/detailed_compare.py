#!/usr/bin/env python3
"""
详细对比RTL和软解输出，分析差异分布
"""
import struct

def read_yuv(filename, width=64, height=64):
    """读取YUV420格式文件"""
    with open(filename, 'rb') as f:
        data = f.read()
    
    y_size = width * height
    uv_width = width // 2
    uv_height = height // 2
    uv_size = uv_width * uv_height
    
    y = list(data[:y_size])
    u = list(data[y_size:y_size+uv_size])
    v = list(data[y_size+uv_size:y_size+2*uv_size])
    
    return y, u, v

def analyze_regions(y_rtl, y_sw, width=64):
    """分析不同区域的差异"""
    # 将64x64分成16个16x16的块
    print("\n16x16块差异分析:")
    print("格式: [块号] RTL平均/软解平均/差异")
    
    for block_row in range(4):
        for block_col in range(4):
            block_id = block_row * 4 + block_col
            
            # 计算该块的像素范围
            start_y = block_row * 16
            end_y = start_y + 16
            start_x = block_col * 16
            end_x = start_x + 16
            
            rtl_sum = 0
            sw_sum = 0
            count = 0
            
            for y in range(start_y, end_y):
                for x in range(start_x, end_x):
                    idx = y * width + x
                    rtl_sum += y_rtl[idx]
                    sw_sum += y_sw[idx]
                    count += 1
            
            rtl_avg = rtl_sum / count
            sw_avg = sw_sum / count
            diff = abs(rtl_avg - sw_avg)
            
            print(f"  [块{block_id:2d}] RTL={rtl_avg:6.2f} SW={sw_avg:6.2f} 差异={diff:6.2f}")

def main():
    # 读取RTL和软解输出
    y_rtl, u_rtl, v_rtl = read_yuv('out/decoded_output.yuv')
    y_sw, u_sw, v_sw = read_yuv('sw_decoder_out_0.yuv')
    
    print("=" * 70)
    print("详细对比分析")
    print("=" * 70)
    
    # Y平面直方图分析
    print("\nY平面直方图 (0-255, 每20一个区间):")
    print("       RTL计数    软解计数    差异")
    
    for i in range(0, 256, 20):
        rtl_count = sum(1 for v in y_rtl if i <= v < i+20)
        sw_count = sum(1 for v in y_sw if i <= v < i+20)
        diff = abs(rtl_count - sw_count)
        print(f"[{i:3d}-{i+19:3d}]  {rtl_count:4d}        {sw_count:4d}        {diff:4d}")
    
    # 区域差异分析
    analyze_regions(y_rtl, y_sw)
    
    # 查找最大的10个差异位置
    print("\nY平面差异最大的10个位置:")
    diffs = []
    for i in range(min(len(y_rtl), len(y_sw))):
        diff = abs(y_rtl[i] - y_sw[i])
        diffs.append((diff, i, y_rtl[i], y_sw[i]))
    
    diffs.sort(reverse=True)
    
    for i, (diff, idx, rtl_val, sw_val) in enumerate(diffs[:10]):
        row = idx // 64
        col = idx % 64
        print(f"  #{i+1}: 位置({row:2d},{col:2d}) RTL={rtl_val:3d} SW={sw_val:3d} 差异={diff:3d}")
    
    # 分析0值分布
    rtl_zeros = sum(1 for v in y_rtl if v == 0)
    sw_zeros = sum(1 for v in y_sw if v == 0)
    print(f"\n0值像素数量:")
    print(f"  RTL: {rtl_zeros} ({rtl_zeros/len(y_rtl)*100:.1f}%)")
    print(f"  软解: {sw_zeros} ({sw_zeros/len(y_sw)*100:.1f}%)")
    
    # 分析高值分布
    rtl_high = sum(1 for v in y_rtl if v > 200)
    sw_high = sum(1 for v in y_sw if v > 200)
    print(f"\n>200值像素数量:")
    print(f"  RTL: {rtl_high} ({rtl_high/len(y_rtl)*100:.1f}%)")
    print(f"  软解: {sw_high} ({sw_high/len(y_sw)*100:.1f}%)")

if __name__ == '__main__':
    main()