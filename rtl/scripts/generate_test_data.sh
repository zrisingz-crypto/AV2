#!/bin/bash
#==============================================================================
# Script Name: generate_test_data.sh
# Purpose: Generate test data files for real decoder simulation
# Usage: ./scripts/generate_test_data.sh
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
echo "Generating Test Data"
echo "========================================"

# Set paths
RTL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="$RTL_DIR/output"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Generate software pixel ROM (gradient pattern 0-255)
echo "Generating software pixel ROM..."
python3 -c "for i in range(4096): print(f'{i%256:02X}')" > "$OUTPUT_DIR/sw_pixel_rom.txt"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Test data generated successfully!${NC}"
    echo "Output: $OUTPUT_DIR/sw_pixel_rom.txt"
    echo "Entries: 4096"
else
    echo -e "${RED}Failed to generate test data!${NC}"
    exit 1
fi