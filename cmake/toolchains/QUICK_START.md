# ARM Cross-Compilation Quick Start Guide

## TL;DR - Most Common Use Cases

### Raspberry Pi 3/4 (32-bit OS)
```bash
# Install compiler
sudo apt-get install gcc-arm-linux-gnueabihf g++-arm-linux-gnueabihf

# Build
./build.sh -T cmake/toolchains/arm-linux-gnueabihf.cmake
```

### Raspberry Pi 3/4/5 (64-bit OS)
```bash
# Install compiler
sudo apt-get install gcc-aarch64-linux-gnu g++-aarch64-linux-gnu

# Build
./build.sh -T cmake/toolchains/aarch64-linux-gnu.cmake
```

### STM32 (Cortex-M4)
```bash
# Install compiler
sudo apt-get install gcc-arm-none-eabi

# Build
./build.sh -T cmake/toolchains/arm-none-eabi.cmake
```

## Command-Line Options

### build.sh options:
- `-T <toolchain_file>` - Specify CMake toolchain file for cross-compilation
- `-N <number>` - Set number of channels (default: 1)
- `-D <level>` - Set logging level (0=none, 1=error, 2=info, 3=debug)

### Examples:
```bash
# Basic cross-compilation
./build.sh -T cmake/toolchains/arm-linux-gnueabihf.cmake

# With 4 channels and debug logging
./build.sh -N 4 -D 3 -T cmake/toolchains/aarch64-linux-gnu.cmake

# STM32 with Cortex-M7
./build.sh -T cmake/toolchains/arm-none-eabi.cmake
# Note: CORTEX_TYPE can be set in the toolchain file or via CMake directly
```

## Testing Your Build

After cross-compiling, the libraries will be in the `build/` directory:
- `build/libpmu_estimator.a` - Static library
- `build/libpmu_estimator.so` - Shared library (not for bare-metal targets)

### Test on Target Device
Copy the library and headers to your target device:
```bash
scp -r build/ pi@raspberrypi.local:~/pmu_estimator/
```

### Test with QEMU (without target device)
```bash
# Install QEMU
sudo apt-get install qemu-user qemu-user-static

# For ARM32
qemu-arm-static -L /usr/arm-linux-gnueabihf/ ./your_test_app

# For ARM64
qemu-aarch64-static -L /usr/aarch64-linux-gnu/ ./your_test_app
```

## Troubleshooting

### "arm-linux-gnueabihf-gcc: not found"
Install the cross-compiler:
```bash
sudo apt-get update
sudo apt-get install gcc-arm-linux-gnueabihf g++-arm-linux-gnueabihf
```

### Build fails with "undefined reference"
For bare-metal targets (STM32), you need to link your final executable with:
- Startup code for your specific MCU
- Linker script for your memory layout
- System libraries (newlib-nano)

### Wrong architecture detected
Verify your toolchain is correctly installed:
```bash
arm-linux-gnueabihf-gcc --version
aarch64-linux-gnu-gcc --version
arm-none-eabi-gcc --version
```

## Next Steps

1. Read the comprehensive documentation: `cmake/toolchains/README.md`
2. Check the main README for general build instructions
3. Review the example programs in `examples/`
4. Test on your target device or using QEMU

## Support Matrix

| Platform | Toolchain File | Status |
|----------|---------------|--------|
| Raspberry Pi 2/3/4 (32-bit) | arm-linux-gnueabihf.cmake | ✅ Tested |
| Raspberry Pi 3/4/5 (64-bit) | aarch64-linux-gnu.cmake | ✅ Tested |
| STM32 (Cortex-M) | arm-none-eabi.cmake | ✅ Tested |
| BeagleBone Black | arm-linux-gnueabihf.cmake | ✅ Compatible |
| NVIDIA Jetson | aarch64-linux-gnu.cmake | ✅ Compatible |
| Nordic nRF | arm-none-eabi.cmake | ✅ Compatible |
