#!/bin/bash
#==============================================================================
# Script Name: run_real_decoder.sh
# Purpose: Run the real AV2 tile decoder simulation
# Usage: ./scripts/run_real_decoder.sh
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
echo "Running Real AV2 Decoder Simulation"
echo "========================================"

# Set paths
RTL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$RTL_DIR"

# Check if simulation executable exists
if [ ! -f "sim_real_decoder" ]; then
    echo -e "${YELLOW}Simulation executable not found. Compiling...${NC}"
    ./scripts/compile_real_decoder.sh
fi

# Run simulation
echo "Starting simulation..."
echo "Output will be saved to sim_real_decoder.log"
vvp tb_ivf_decode_real 2>&1 | tee sim_real_decoder.log

if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo -e "${GREEN}Simulation completed successfully!${NC}"
    echo "Log: sim_real_decoder.log"
else
    echo -e "${YELLOW}Simulation finished with exit code ${PIPESTATUS[0]}${NC}"
    echo "Check sim_real_decoder.log for details"
fi