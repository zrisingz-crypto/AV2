#!/bin/bash
#==============================================================================
# Script Name: compile_real_decoder.sh
# Purpose: Compile the real AV2 tile decoder with testbench
# Usage: ./scripts/compile_real_decoder.sh
# Author: AV2 RTL Team
# Date: 2026-02-12
#==============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "========================================"
echo "Compiling Real AV2 Decoder"
echo "========================================"

# Set paths
RTL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$RTL_DIR"

# Compile with iverilog
echo "Compiling RTL modules..."
iverilog -g2012 \
    -o sim_real_decoder \
    -I. \
    tb_tile_decoder_real.v \
    av2_tile_decoder_real.v \
    av2_entropy_decoder_real.v \
    av2_coeff_decoder_real.v \
    av2_intra_prediction_real_fixed.v \
    av2_inverse_transform_real_fixed.v \
    2>&1 | tee compile_real_decoder.log

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Compilation successful!${NC}"
    echo "Output: sim_real_decoder"
    echo "Log: compile_real_decoder.log"
else
    echo -e "${RED}Compilation failed!${NC}"
    exit 1
fi