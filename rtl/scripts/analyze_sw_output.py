#!/usr/bin/env python3

# 读取软解码器输出
with open('sw_decoder_out_0.yuv', 'rb') as f:
    data = f.read()

# 提取各个平面
y_plane = [list(data[i*64:(i+1)*64]) for i in range(64)]
u_plane = [list(data[4096+i*32:4096+(i+1)*32]) for i in range(32)]
v_plane = [list(data[5120+i*32:5120+(i+1)*32]) for i in range(32)]

print("=" * 60)
print("软件解码器输出分析")
print("=" * 60)

print("\n[Y平面统计]")
y_flat = [pixel for row in y_plane for pixel in row]
print(f"  最小值: {min(y_flat)}")
print(f"  最大值: {max(y_flat)}")
print(f"  平均值: {sum(y_flat)/len(y_flat):.2f}")

# 检查左上角区域
print(f"\n[Y平面左上角16x16区域值]")
for i in range(16):
    if i % 4 == 0:
        print(f"  行{i:2d}: {y_plane[i][:16]}")
    else:
        print(f"        {y_plane[i][:16]}")

# 检查是否有明显的模式
print(f"\n[U平面前8x8区域]")
for i in range(8):
    print(f"  行{i}: {u_plane[i][:8]}")

print(f"\n[V平面前8x8区域]")
for i in range(8):
    print(f"  行{i}: {v_plane[i][:8]}")

# 检查边缘特征
print(f"\n[边缘特征分析]")
print(f"  Y平面左上角4x4:")
for i in range(4):
    print(f"    {y_plane[i][:4]}")
print(f"  Y平面右下角4x4:")
for i in range(60, 64):
    print(f"    {y_plane[i][60:64]}")
