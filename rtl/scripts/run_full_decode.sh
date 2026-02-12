#!/bin/bash
# AV2 RTL 硬件解码完整流程一键执行脚本
# 用法: ./run_full_decode.sh

set -e

echo "============================================"
echo "AV2 RTL 硬件解码完整流程"
echo "============================================"

# 检查依赖
echo ""
echo "[0/5] 检查依赖..."
if ! command -v iverilog &> /dev/null; then
    echo "错误: iverilog 未安装!"
    echo "请使用 'brew install icarus-verilog' 安装"
    exit 1
fi

if ! command -v vvp &> /dev/null; then
    echo "错误: vvp 未安装!"
    exit 1
fi

if ! python3 -c "from PIL import Image" 2>/dev/null; then
    echo "警告: PIL/Pillow 未安装，将尝试安装..."
    pip3 install Pillow
fi

echo "依赖检查通过!"

# 步骤 1: 编译
echo ""
echo "[1/5] 编译 RTL 代码..."
iverilog -g2012 -o tb_tile_decoder_real \
    tb_tile_decoder_real.v \
    av2_tile_decoder_real.v \
    av2_intra_prediction_real_fixed.v \
    av2_inverse_transform_real_fixed.v \
    av2_coeff_decoder_real.v \
    av2_entropy_decoder_real.v 2>&1 | tee compile.log

if grep -q "error" compile.log; then
    echo "错误: 编译失败!"
    exit 1
fi

echo "编译完成!"

# 步骤 2: 运行仿真
echo ""
echo "[2/5] 运行仿真..."
rm -f rtl_real_output.txt
vvp tb_tile_decoder_real -l sim_real.log 2>&1 | tail -20

# 步骤 3: 检查输出
echo ""
echo "[3/5] 检查输出..."
if [ ! -f rtl_real_output.txt ]; then
    echo "错误: rtl_real_output.txt 未生成!"
    exit 1
fi

LINE_COUNT=$(wc -l < rtl_real_output.txt)
echo "RTL 输出文件: rtl_real_output.txt ($LINE_COUNT 行)"

# 检查是否有 'x' 值
X_COUNT=$(grep -c "= x" rtl_real_output.txt || true)
if [ "$X_COUNT" -gt 0 ]; then
    echo "警告: 发现 $X_COUNT 个未知值 (x)"
else
    echo "状态: 所有像素值有效 (无 'x' 值)"
fi

# 步骤 4: 转 YUV
echo ""
echo "[4/5] 转换为 YUV..."
python3 convert_rtl_to_yuv.py rtl_real_output.txt hw_decoder_out_real.yuv

YUV_SIZE=$(wc -c < hw_decoder_out_real.yuv)
echo "YUV 文件: hw_decoder_out_real.yuv ($YUV_SIZE 字节)"

# 验证 YUV 大小 (64x64 YUV420 = 6144 字节)
if [ "$YUV_SIZE" -ne 6144 ]; then
    echo "警告: YUV 文件大小异常，期望 6144 字节"
fi

# 步骤 5: 转 PNG
echo ""
echo "[5/5] 转换为 PNG..."
python3 yuv_to_png.py hw_decoder_out_real.yuv hw_decoder_out_real.png
echo "PNG 文件: hw_decoder_out_real.png"

# 与软件解码对比（如果存在）
echo ""
echo "[额外] 与软件解码对比..."
if [ -f output/sw_output.yuv ]; then
    python3 output/compare_real_outputs.py rtl_real_output.txt output/sw_output.yuv
else
    echo "软件解码输出不存在，跳过对比"
fi

# 完成
echo ""
echo "============================================"
echo "全部完成!"
echo "============================================"
echo "输出文件:"
echo "  - RTL 文本: rtl_real_output.txt ($(wc -l < rtl_real_output.txt) 行)"
echo "  - YUV 文件: hw_decoder_out_real.yuv ($(wc -c < hw_decoder_out_real.yuv) 字节)"
echo "  - PNG 图像: hw_decoder_out_real.png"
echo "  - 仿真日志: sim_real.log"
echo "============================================"
echo ""
echo "查看图像:"
echo "  open hw_decoder_out_real.png"
