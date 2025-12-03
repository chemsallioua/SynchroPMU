# ARM Cross-Compilation Toolchain Files

This directory contains CMake toolchain files for cross-compiling the PmuEstimator library to various ARM architectures.

## Available Toolchains

### 1. arm-linux-gnueabihf.cmake
**Target**: 32-bit ARM Linux with hardware floating-point
**Suitable for**:
- Raspberry Pi 2/3/4 (running 32-bit Raspberry Pi OS)
- BeagleBone Black
- Other ARM Cortex-A devices running Linux

**Prerequisites**:
```bash
sudo apt-get install gcc-arm-linux-gnueabihf g++-arm-linux-gnueabihf
```

**Usage**:
```bash
cmake -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-linux-gnueabihf.cmake ..
```

### 2. aarch64-linux-gnu.cmake
**Target**: 64-bit ARM Linux (AArch64)
**Suitable for**:
- Raspberry Pi 3/4/5 (running 64-bit Raspberry Pi OS)
- NVIDIA Jetson series
- Other ARM64 devices running Linux

**Prerequisites**:
```bash
sudo apt-get install gcc-aarch64-linux-gnu g++-aarch64-linux-gnu
```

**Usage**:
```bash
cmake -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/aarch64-linux-gnu.cmake ..
```

### 3. arm-none-eabi.cmake
**Target**: Bare-metal ARM Cortex-M microcontrollers
**Suitable for**:
- STM32 series (STM32F4, STM32H7, etc.)
- Nordic nRF series
- Other ARM Cortex-M based microcontrollers

**Prerequisites**:
```bash
sudo apt-get install gcc-arm-none-eabi
```

**Usage**:
```bash
# For Cortex-M4 (default)
cmake -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-none-eabi.cmake ..

# For Cortex-M7
cmake -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-none-eabi.cmake \
      -DCORTEX_TYPE=M7 ..

# For Cortex-M3
cmake -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-none-eabi.cmake \
      -DCORTEX_TYPE=M3 ..
```

Supported CORTEX_TYPE values: M0, M0PLUS, M3, M4, M7, M33

## Complete Build Examples

### Example 1: Build for Raspberry Pi 3 (32-bit)
```bash
# Install cross-compiler
sudo apt-get install gcc-arm-linux-gnueabihf g++-arm-linux-gnueabihf

# Build using build.sh
./build.sh -T cmake/toolchains/arm-linux-gnueabihf.cmake

# Or using CMake directly
mkdir build && cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/arm-linux-gnueabihf.cmake ..
make PmuEstimatorStatic
make PmuEstimatorShared
```

### Example 2: Build for Raspberry Pi 4 (64-bit)
```bash
# Install cross-compiler
sudo apt-get install gcc-aarch64-linux-gnu g++-aarch64-linux-gnu

# Build with custom settings
./build.sh -T cmake/toolchains/aarch64-linux-gnu.cmake -N 4 -D 2

# Or using CMake directly
mkdir build && cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/aarch64-linux-gnu.cmake \
      -DNUM_CHANLS=4 \
      -DLOGGING_LEVEL=2 ..
make
```

### Example 3: Build for STM32F4 (Cortex-M4)
```bash
# Install cross-compiler
sudo apt-get install gcc-arm-none-eabi

# Build for Cortex-M4
./build.sh -T cmake/toolchains/arm-none-eabi.cmake

# Or using CMake directly
mkdir build && cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/arm-none-eabi.cmake \
      -DCORTEX_TYPE=M4 ..
make PmuEstimatorStatic  # Note: Shared libraries not supported for bare-metal
```

## Custom Toolchain Files

You can create your own toolchain file for specific ARM platforms. Use the provided files as templates and adjust:
- Compiler paths
- CPU-specific flags
- FPU configuration
- Target system directories

## Troubleshooting

### Cross-compiler not found
Make sure the cross-compiler is installed and in your PATH:
```bash
which arm-linux-gnueabihf-gcc
which aarch64-linux-gnu-gcc
which arm-none-eabi-gcc
```

### Linking errors on embedded targets
For bare-metal ARM targets, you may need to provide additional linker scripts and startup code specific to your microcontroller. The toolchain file provides a generic configuration that works for building the library, but final linking into an executable may require platform-specific setup.

### Testing cross-compiled binaries
You can test Linux ARM binaries using QEMU:
```bash
# Install QEMU
sudo apt-get install qemu-user qemu-user-static

# Run ARM binary
qemu-arm-static -L /usr/arm-linux-gnueabihf/ ./your_arm_binary

# Run ARM64 binary
qemu-aarch64-static -L /usr/aarch64-linux-gnu/ ./your_arm64_binary
```

## Additional Resources

- [CMake Cross Compiling](https://cmake.org/cmake/help/latest/manual/cmake-toolchains.7.html#cross-compiling)
- [ARM GCC Toolchains](https://developer.arm.com/tools-and-software/open-source-software/developer-tools/gnu-toolchain)
- [Raspberry Pi Cross-Compilation](https://www.raspberrypi.com/documentation/computers/linux_kernel.html#cross-compiling-the-kernel)
