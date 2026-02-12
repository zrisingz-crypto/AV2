#!/bin/bash
#==============================================================================
# Script Name: compile_ivf_real_decoder.sh
# Purpose: Compile the real IVF decoder testbench
# Usage: ./scripts/compile_ivf_real_decoder.sh
# Author: AV2 RTL Team
# Date: 2026-02-12
#==============================================================================

set -e  # Exit on error

echo "========================================"
echo "Compiling Real IVF Decoder"
echo "========================================"

# Change to rtl directory
cd "$(dirname "$0")/.."

# Clean previous build
echo "[INFO] Cleaning previous build..."
rm -f tb_ivf_decode_real

# Compile with Icarus Verilog
echo "[INFO] Compiling testbench..."
iverilog -o tb_ivf_decode_real \
    -g2012 \
    -Wall \
    tb_ivf_decode_real.v \
    ivf_bitstream_data.v \
    av2_tile_decoder_real.v \
    av2_entropy_decoder_real.v \
    av2_coeff_decoder_real.v \
    av2_intra_prediction_real_fixed.v \
    av2_inverse_transform_real_fixed.v \
    av2_context_model.v \
    av2_frame_buffer_ctrl.v \
    av2_output_ctrl.v

# Check compilation result
if [ $? -eq 0 ]; then
    echo "========================================"
    echo "✓ Compilation successful!"
    echo "========================================"
    echo "[INFO] Executable: tb_ivf_decode_real"
    echo "[INFO] Run with: vvp tb_ivf_decode_real"
else
    echo "========================================"
    echo "✗ Compilation failed!"
    echo "========================================"
    exit 1
fi