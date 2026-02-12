#!/bin/bash

#==============================================================================
# Run Unit Tests for RTL Modules
#==============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Directories
RTL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="$RTL_DIR/output"
TEST_RESULTS_DIR="$OUTPUT_DIR/test_results"

# Create directories
mkdir -p "$TEST_RESULTS_DIR"

echo "========================================="
echo "RTL Unit Tests Runner"
echo "========================================="
echo "RTL Directory: $RTL_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo ""

# Test results summary
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to compile and run a test
run_test() {
    local test_name=$1
    local testbench=$2
    local dut=$3
    local test_result_file="$TEST_RESULTS_DIR/${test_name}_result.txt"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    echo -e "${YELLOW}Running Test: $test_name${NC}"
    echo "  Testbench: $testbench"
    echo "  DUT: $dut"
    
    # Compile
    echo "  Compiling..."
    local compile_cmd="iverilog -g2012 -o \"$OUTPUT_DIR/${test_name}_sim\" \
        \"$RTL_DIR/$testbench\" \
        \"$RTL_DIR/$dut\""
    
    if ! eval "$compile_cmd" 2>&1 | tee "$test_result_file"; then
        echo -e "${RED}✗ FAILED: Compilation error${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Run simulation
    echo "  Running simulation..."
    if ! vvp "$OUTPUT_DIR/${test_name}_sim" 2>&1 | tee -a "$test_result_file"; then
        echo -e "${RED}✗ FAILED: Simulation error${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Check if test passed
    if grep -q "ALL TESTS PASSED" "$test_result_file"; then
        echo -e "${GREEN}✓ PASSED${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    else
        echo -e "${RED}✗ FAILED: Tests did not pass${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
}

# Run all tests
echo ""
echo "========================================="
echo "Running Tests"
echo "========================================="
echo ""

# Test 1: Motion Compensation
run_test "motion_compensation" \
    "tb_motion_compensation.v" \
    "av2_motion_compensation_real.v" || true

echo ""

# Test 2: MV Decoder
run_test "mv_decoder" \
    "tb_mv_decoder.v" \
    "av2_mv_decoder_real.v" || true

echo ""

# Test 3: Deblocking Filter
run_test "deblocking_filter" \
    "tb_deblocking_filter.v" \
    "av2_deblocking_filter_real.v" || true

echo ""

# Test 4: CDEF Filter
run_test "cdef_filter" \
    "tb_cdef_filter.v" \
    "av2_cdef_filter_real.v" || true

echo ""
echo "========================================="
echo "Test Summary"
echo "========================================="
echo "Total Tests: $TOTAL_TESTS"
echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"
echo -e "${RED}Failed: $FAILED_TESTS${NC}"
echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}✓ ALL TESTS PASSED${NC}"
    exit 0
else
    echo -e "${RED}✗ SOME TESTS FAILED${NC}"
    echo ""
    echo "Check test results in: $TEST_RESULTS_DIR"
    exit 1
fi