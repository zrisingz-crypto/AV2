#!/usr/bin/env python3

# 读取RTL解码器输出
with open('out/decoded_output.yuv', 'rb') as f:
    rtl_data = f.read()

# 读取软解码器输出
with open('sw_decoder_out_0.yuv', 'rb') as f:
    sw_data = f.read()

print("=" * 60)
print("RTL vs 软解码器输出对比")
print("=" * 60)

print(f"\n[文件大小]")
print(f"  RTL输出:  {len(rtl_data)} 字节")
print(f"  软解输出:  {len(sw_data)} 字节")

# 提取Y平面
if len(rtl_data) >= 4096 and len(sw_data) >= 4096:
    rtl_y = rtl_data[:4096]
    sw_y = sw_data[:4096]
    
    print(f"\n[Y平面统计]")
    
    rtl_y_flat = list(rtl_y)
    sw_y_flat = list(sw_y)
    
    print(f"  RTL: 最小={min(rtl_y_flat)}, 最大={max(rtl_y_flat)}, 平均={sum(rtl_y_flat)/len(rtl_y_flat):.2f}")
    print(f"  软解: 最小={min(sw_y_flat)}, 最大={max(sw_y_flat)}, 平均={sum(sw_y_flat)/len(sw_y_flat):.2f}")
    
    # 计算差异
    diff_count = 0
    diff_sum = 0
    diff_max = 0
    
    for i in range(4096):
        diff = abs(rtl_y[i] - sw_y[i])
        if diff > 0:
            diff_count += 1
            diff_sum += diff
            if diff > diff_max:
                diff_max = diff
    
    print(f"\n[Y平面差异分析]")
    print(f"  不同像素数: {diff_count} / 4096 ({diff_count/4096*100:.2f}%)")
    if diff_count > 0:
        print(f"  平均差异: {diff_sum/diff_count:.2f}")
        print(f"  最大差异: {diff_max}")
    
    # 检查关键区域
    print(f"\n[左上角4x4区域对比]")
    print("  RTL:")
    for row in range(4):
        print(f"    {list(rtl_y[row*64:row*64+4])}")
    print("  软解:")
    for row in range(4):
        print(f"    {list(sw_y[row*64:row*64+4])}")
    
    print(f"\n[右下角4x4区域对比]")
    print("  RTL:")
    for row in range(60, 64):
        print(f"    {list(rtl_y[row*64+60:row*64+64])}")
    print("  软解:")
    for row in range(60, 64):
        print(f"    {list(sw_y[row*64+60:row*64+64])}")

# 检查UV平面
if len(rtl_data) >= 6144 and len(sw_data) >= 6144:
    rtl_u = rtl_data[4096:5120]
    rtl_v = rtl_data[5120:6144]
    sw_u = sw_data[4096:5120]
    sw_v = sw_data[5120:6144]
    
    print(f"\n[U平面统计]")
    print(f"  RTL: 最小={min(rtl_u)}, 最大={max(rtl_u)}, 平均={sum(rtl_u)/len(rtl_u):.2f}")
    print(f"  软解: 最小={min(sw_u)}, 最大={max(sw_u)}, 平均={sum(sw_u)/len(sw_u):.2f}")
    
    print(f"\n[V平面统计]")
    print(f"  RTL: 最小={min(rtl_v)}, 最大={max(rtl_v)}, 平均={sum(rtl_v)/len(rtl_v):.2f}")
    print(f"  软解: 最小={min(sw_v)}, 最大={max(sw_v)}, 平均={sum(sw_v)/len(sw_v):.2f}")
    
    # UV平面差异
    u_diff = sum(1 for i in range(1024) if rtl_u[i] != sw_u[i])
    v_diff = sum(1 for i in range(1024) if rtl_v[i] != sw_v[i])
    
    print(f"\n[UV平面差异]")
    print(f"  U平面不同像素: {u_diff} / 1024 ({u_diff/1024*100:.2f}%)")
    print(f"  V平面不同像素: {v_diff} / 1024 ({v_diff/1024*100:.2f}%)")

print(f"\n[整体评估]")
if len(rtl_data) >= 6144:
    total_diff = 0
    for i in range(min(len(rtl_data), len(sw_data))):
        if rtl_data[i] != sw_data[i]:
            total_diff += 1
    match_rate = (min(len(rtl_data), len(sw_data)) - total_diff) / min(len(rtl_data), len(sw_data)) * 100
    print(f"  整体匹配率: {match_rate:.2f}%")
    print(f"  不同字节数: {total_diff}")