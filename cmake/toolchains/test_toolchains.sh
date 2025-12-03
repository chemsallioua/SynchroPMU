#!/bin/bash
# Test script to verify toolchain files are correctly configured
# Note: This requires cross-compilers to be installed
#
# Usage: ./test_toolchains.sh [toolchain_name]
#   Where toolchain_name is one of: arm32, arm64, stm32, all

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "Repository root: $REPO_ROOT"
echo "Toolchains directory: $SCRIPT_DIR"

test_toolchain() {
    local name=$1
    local toolchain_file=$2
    local extra_flags=$3
    
    echo ""
    echo "========================================"
    echo "Testing: $name"
    echo "Toolchain: $toolchain_file"
    echo "========================================"
    
    # Create test build directory
    local build_dir="$REPO_ROOT/cmake_test_$name"
    rm -rf "$build_dir"
    mkdir -p "$build_dir"
    cd "$build_dir"
    
    # Try to configure with the toolchain
    if cmake -DCMAKE_TOOLCHAIN_FILE="$toolchain_file" $extra_flags ..; then
        echo "✓ Configuration successful for $name"
        
        # Try to build
        if make PmuEstimatorStatic; then
            echo "✓ Build successful for $name"
            echo "✓ Library created: $(ls -lh libpmu_estimator.a 2>/dev/null || echo 'not found')"
        else
            echo "✗ Build failed for $name"
            return 1
        fi
    else
        echo "✗ Configuration failed for $name"
        echo "  (This is expected if the cross-compiler is not installed)"
        return 1
    fi
    
    # Clean up
    cd "$REPO_ROOT"
    rm -rf "$build_dir"
}

test_arm32() {
    test_toolchain "ARM32-Linux" \
                   "$SCRIPT_DIR/arm-linux-gnueabihf.cmake" \
                   ""
}

test_arm64() {
    test_toolchain "ARM64-Linux" \
                   "$SCRIPT_DIR/aarch64-linux-gnu.cmake" \
                   ""
}

test_stm32() {
    test_toolchain "STM32-M4" \
                   "$SCRIPT_DIR/arm-none-eabi.cmake" \
                   "-DCORTEX_TYPE=M4"
}

# Check what test to run
TEST_TARGET="${1:-all}"

case "$TEST_TARGET" in
    arm32)
        test_arm32
        ;;
    arm64)
        test_arm64
        ;;
    stm32)
        test_stm32
        ;;
    all)
        echo "Testing all toolchains..."
        failed=0
        
        test_arm32 || ((failed++))
        test_arm64 || ((failed++))
        test_stm32 || ((failed++))
        
        echo ""
        echo "========================================"
        echo "Summary"
        echo "========================================"
        if [ $failed -eq 0 ]; then
            echo "✓ All tests passed!"
        else
            echo "✗ $failed test(s) failed"
            echo "  (This is expected if cross-compilers are not installed)"
        fi
        ;;
    *)
        echo "Usage: $0 [arm32|arm64|stm32|all]"
        exit 1
        ;;
esac
