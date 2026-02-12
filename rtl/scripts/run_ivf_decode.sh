#!/bin/bash

#==============================================================================
# IVF Bitstream Decode Simulation Script
# Decodes real IVF file and outputs YUV
#==============================================================================

echo "========================================"
echo "IVF Bitstream Decode Simulation"
echo "========================================"

# 检查 Icarus Verilog 是否安装
if ! command -v iverilog &> /dev/null; then
    echo "Error: Icarus Verilog not found!"
    echo "Please install it using:"
    echo "  macOS: brew install icarus-verilog"
    echo "  Linux: sudo apt-get install iverilog"
    exit 1
fi

# 进入 RTL 目录
cd "$(dirname "$0")"

# 创建输出目录
mkdir -p out
mkdir -p log

# 检查IVF bitstream数据文件是否存在
if [ ! -f "ivf_bitstream_data.v" ]; then
    echo "Error: ivf_bitstream_data.v not found!"
    echo "Please run: python3 ivf_to_verilog.py"
    exit 1
fi

# 检查hex文件是否存在
if [ ! -f "ivf_bitstream_data_frame_0.hex" ] || [ ! -f "ivf_bitstream_data_frame_1.hex" ]; then
    echo "Error: hex data files not found!"
    echo "Please run: python3 ivf_to_verilog.py"
    exit 1
fi

echo ""
echo "[1/3] Compiling Verilog files..."

# 编译所有 Verilog 文件
iverilog -o ivf_decode_sim \
    -g2012 \
    -Wall \
    av2_decoder_top.v \
    av2_tile_decoder_real.v \
    av2_entropy_decoder_real.v \
    av2_coeff_decoder_real.v \
    av2_inverse_transform_real_fixed.v \
    av2_intra_prediction_real_fixed.v \
    av2_context_model.v \
    av2_mv_decoder.v \
    av2_motion_compensation.v \
    av2_deblocking_filter.v \
    av2_cdef_filter.v \
    av2_obu_parser.v \
    av2_frame_header_parser.v \
    av2_frame_buffer_ctrl.v \
    av2_output_ctrl.v \
    tb_ivf_decode.v

if [ $? -ne 0 ]; then
    echo "✗ Compilation failed!"
    exit 1
fi

echo "✓ Compilation successful!"

echo ""
echo "[2/3] Running simulation..."

# 运行仿真，重定向日志
vvp ivf_decode_sim 2>&1 | tee log/simulation.log

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "✗ Simulation failed!"
    exit 1
fi

echo "✓ Simulation completed!"

echo ""
echo "[3/3] Checking results..."

# 检查是否生成了YUV输出文件
if [ -f "out/decoded_output.yuv" ]; then
    YUV_SIZE=$(stat -f%z "out/decoded_output.yuv" 2>/dev/null || stat -c%s "out/decoded_output.yuv" 2>/dev/null)
    echo "✓ YUV output generated: out/decoded_output.yuv"
    echo "  File size: $YUV_SIZE bytes"
    
    # 计算64x64 YUV420应该的大小
    # Y = 64*64 = 4096 bytes
    # U = 32*32 = 1024 bytes
    # V = 32*32 = 1024 bytes
    # Total per frame = 6144 bytes
    # 2 frames = 12288 bytes
    EXPECTED_SIZE=12288
    
    if [ "$YUV_SIZE" -eq "$EXPECTED_SIZE" ]; then
        echo "  ✓ Size matches expected ($EXPECTED_SIZE bytes)"
    else
        echo "  ⚠ Size differs from expected ($EXPECTED_SIZE bytes)"
    fi
else
    echo "⚠ YUV output file not found"
fi

echo ""
echo "========================================"
echo "Simulation Complete!"
echo "========================================"
echo ""
echo "Output files:"
echo "  YUV: out/decoded_output.yuv"
echo "  Log: log/simulation.log"
echo ""
echo "To view the decoded YUV video:"
echo "  ffplay -f rawvideo -video_size 64x64 -pixel_format yuv420p out/decoded_output.yuv"
echo ""
echo "To convert to a viewable format:"
echo "  ffmpeg -f rawvideo -video_size 64x64 -pixel_format yuv420p -i out/decoded_output.yuv -pix_fmt rgb24 -f image2 out/decoded_frame_%d.png"
