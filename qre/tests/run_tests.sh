#!/bin/bash
# Run all tests for Trotter error coefficient computation

echo "========================================================================"
echo "TROTTER ERROR COEFFICIENT TESTS"
echo "========================================================================"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

FAILED=0
PASSED=0

run_test() {
    local test_file=$1
    local test_name=$2

    echo "----------------------------------------"
    echo "Running: $test_name"
    echo "----------------------------------------"

    if python3.11 "$test_file" > /tmp/test_output_$$.txt 2>&1; then
        echo "✅ PASSED"
        ((PASSED++))
    else
        echo "❌ FAILED"
        echo "Output:"
        cat /tmp/test_output_$$.txt
        ((FAILED++))
    fi
    echo ""
    rm -f /tmp/test_output_$$.txt
}

# Run all test suites
run_test "test_exact_computation.py" "Exact Computation (analytical + edge cases)"
run_test "test_modes.py" "Mode Selection and Integration"

echo "========================================================================"
echo "TEST SUMMARY"
echo "========================================================================"
echo "  Passed: $PASSED"
echo "  Failed: $FAILED"
echo "========================================================================"

if [ $FAILED -eq 0 ]; then
    echo "✅ ALL TESTS PASSED!"
    exit 0
else
    echo "❌ SOME TESTS FAILED"
    exit 1
fi
