#!/usr/bin/env python3
"""
分析软件解码器的输出模式，理解其预测方式
"""
import struct

# 分析软件解码器输出的特征
with open('sw_decoder_out_0.yuv', 'rb') as f:
    y_data = struct.unpack('4096B', f.read(4096))
    u_data = struct.unpack('1024B', f.read(1024))
    v_data = struct.unpack('1024B', f.read(1024))

# 检查是否有重复模式
print('Y平面模式分析:')
unique_values = len(set(y_data))
print(f'  唯一值数量: {unique_values}/4096')

# 检查水平梯度
horizontal_grad = []
for r in range(64):
    for c in range(63):
        horizontal_grad.append(abs(int(y_data[r*64 + c]) - int(y_data[r*64 + c + 1])))
print(f'  水平梯度平均值: {sum(horizontal_grad)/len(horizontal_grad):.2f}')
print(f'  水平梯度最大值: {max(horizontal_grad)}')

# 检查垂直梯度
vertical_grad = []
for r in range(63):
    for c in range(64):
        vertical_grad.append(abs(int(y_data[r*64 + c]) - int(y_data[(r+1)*64 + c])))
print(f'  垂直梯度平均值: {sum(vertical_grad)/len(vertical_grad):.2f}')
print(f'  垂直梯度最大值: {max(vertical_grad)}')

# 检查是否是DC预测（值应该相对平滑）
if unique_values < 100:
    print('  可能是DC预测或简单模式')
elif max(horizontal_grad) < 10 and max(vertical_grad) < 10:
    print('  可能是平滑预测')
else:
    print('  可能是真实视频内容或复杂预测模式')

# 检查块内平滑度
print('\n16x16块内平滑度分析:')
for block_row in range(4):
    for block_col in range(4):
        block_grads = []
        for r in range(16):
            row = block_row * 16 + r
            for c in range(15):
                col = block_col * 16 + c
                block_grads.append(abs(int(y_data[row*64 + col]) - int(y_data[row*64 + col + 1])))
        avg_grad = sum(block_grads) / len(block_grads)
        print(f'  Block({block_row},{block_col}): 平均梯度={avg_grad:.2f}')

# 尝试检测是否是真实的帧内预测（如DC、PLANAR、ANGULAR等）
print('\n预测模式推测:')

# DC预测检测：块内值应该接近平均值
dc_like_blocks = 0
for block_row in range(4):
    for block_col in range(4):
        block_values = []
        for r in range(16):
            for c in range(16):
                idx = (block_row*16 + r) * 64 + (block_col*16 + c)
                block_values.append(y_data[idx])
        avg = sum(block_values) / len(block_values)
        variance = sum((x - avg)**2 for x in block_values) / len(block_values)
        if variance < 100:  # 低方差，可能是DC预测
            dc_like_blocks += 1
            print(f'  Block({block_row},{block_col}): 可能是DC预测 (方差={variance:.2f})')

print(f'\n总共 {dc_like_blocks}/16 个块可能是DC预测')

# 检查是否是PLANAR预测（平滑过渡）
print('\n检查是否有PLANAR预测特征:')
# PLANAR预测在块内应该有平滑的水平/垂直渐变
planar_like = 0
for block_row in range(4):
    for block_col in range(4):
        # 检查水平方向是否有平滑变化
        row_diffs = []
        for r in range(8, 16):  # 检查下半部分
            row = block_row * 16 + r
            left_avg = sum(y_data[row*64 + block_col*16 + c] for c in range(8)) / 8
            right_avg = sum(y_data[row*64 + block_col*16 + c] for c in range(8, 16)) / 8
            row_diffs.append(abs(left_avg - right_avg))
        
        if len(row_diffs) > 0 and sum(row_diffs) / len(row_diffs) < 20:
            planar_like += 1

print(f'  {planar_like}/16 个块可能有PLANAR预测特征')

# 最终结论
print('\n' + '='*60)
print('结论:')
if dc_like_blocks >= 12:
    print('  大多数块使用DC预测')
elif planar_like >= 8:
    print('  很多块使用PLANAR预测')
else:
    print('  使用混合预测模式或真实视频内容')
print('='*60)