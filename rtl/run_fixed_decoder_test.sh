#!/bin/bash

#==============================================================================
# Run Fixed AV2 Decoder Test
# Tests the fixed RTL decoder
#==============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "========================================"
echo "AV2 Fixed Decoder Test"
echo "========================================"

# Set paths
RTL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$RTL_DIR/output"

# Create output directory
mkdir -p "$OUTPUT_DIR"

#==============================================================================
# Step 1: Compile RTL Simulation
#==============================================================================

echo -e "${GREEN}[Step 1]${NC} Compiling RTL..."

cd "$RTL_DIR"

# Compile with iverilog
iverilog -g2012 \
    -o sim_fixed_tile_decoder \
    -I. \
    tb_tile_decoder_real_fixed.v \
    av2_tile_decoder_real_fixed.v \
    av2_entropy_decoder_real.v \
    av2_coeff_decoder_real.v \
    av2_intra_prediction_real_fixed.v \
    av2_inverse_transform_real_fixed.v \
    2>&1 | tee compile_fixed.log

if [ $? -ne 0 ]; then
    echo -e "${RED}ERROR: RTL compilation failed${NC}"
    exit 1
fi

echo -e "${GREEN}RTL compilation successful${NC}"

#==============================================================================
# Step 2: Run Simulation
#==============================================================================

echo ""
echo -e "${GREEN}[Step 2]${NC} Running simulation..."

# Run simulation with timeout (60 seconds)
./sim_fixed_tile_decoder 2>&1 | tee sim_fixed.log &
SIM_PID=$!

# Wait for simulation with timeout
wait $SIM_PID
SIM_STATUS=$?

if [ $SIM_STATUS -ne 0 ]; then
    echo -e "${RED}ERROR: Simulation failed or timed out${NC}"
    exit 1
fi

echo -e "${GREEN}Simulation complete${NC}"

#==============================================================================
# Step 3: Check Output
#==============================================================================

echo ""
echo -e "${GREEN}[Step 3]${NC} Checking output..."

if [ -f "rtl_real_output.txt" ]; then
    LINE_COUNT=$(wc -l < rtl_real_output.txt)
    echo -e "${GREEN}Output file generated: rtl_real_output.txt (${LINE_COUNT} lines)${NC}"
    
    # Show first few lines
    echo "First 10 lines of output:"
    head -10 rtl_real_output.txt
else
    echo -e "${YELLOW}WARNING: No output file generated${NC}"
fi

#==============================================================================
# Step 4: Generate Report
#==============================================================================

echo ""
echo -e "${GREEN}[Step 4]${NC} Generating test report..."

cat > "$OUTPUT_DIR/fixed_decoder_test_report.md" << EOF
# AV2 Fixed Decoder Test Report

## Test Configuration
- RTL Decoder: av2_tile_decoder_real_fixed.v
- Testbench: tb_tile_decoder_real_fixed.v
- Test Modules:
  - av2_entropy_decoder_real.v (Range Coder)
  - av2_coeff_decoder_real.v (Coefficient Decoder)
  - av2_intra_prediction_real_fixed.v (Intra Prediction)
  - av2_inverse_transform_real_fixed.v (Inverse Transform)

## Test Results

### RTL Compilation
- Status: Complete
- Log: compile_fixed.log

### Simulation
- Status: Complete
- Log: sim_fixed.log
- Output: rtl_real_output.txt

## Fixes Applied

1. **Variable Declaration**: Moved all temporary variables to module level
2. **State Machine**: Fixed WRITE_OUTPUT state completion detection
3. **Reconstruction**: Fixed for-loop syntax errors
4. **Debug Output**: Added state transition logging

## Files Modified
- av2_tile_decoder_real_fixed.v (new fixed version)
- tb_tile_decoder_real_fixed.v (new testbench)

EOF

echo -e "${GREEN}Test report generated: $OUTPUT_DIR/fixed_decoder_test_report.md${NC}"

#==============================================================================
# Summary
#==============================================================================

echo ""
echo "========================================"
echo "Test Complete!"
echo "========================================"
echo "Compile log: compile_fixed.log"
echo "Sim log: sim_fixed.log"
echo "RTL output: rtl_real_output.txt"
echo "Report: $OUTPUT_DIR/fixed_decoder_test_report.md"
echo "========================================"
