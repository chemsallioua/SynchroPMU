# CMake Toolchain file for cross-compiling to ARM Cortex-M (bare-metal)
# Suitable for STM32, Nordic nRF, and other ARM Cortex-M microcontrollers
#
# Usage:
#   cmake -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-none-eabi.cmake \
#         -DCORTEX_TYPE=M4 ..
#
# Set CORTEX_TYPE to: M0, M0PLUS, M3, M4, M7, M33, etc.
#

set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

# Specify the cross compiler
set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_ASM_COMPILER arm-none-eabi-gcc)
set(CMAKE_AR arm-none-eabi-ar)
set(CMAKE_OBJCOPY arm-none-eabi-objcopy)
set(CMAKE_OBJDUMP arm-none-eabi-objdump)
set(CMAKE_SIZE arm-none-eabi-size)

# Skip compiler test (bare-metal doesn't link without startup code)
set(CMAKE_C_COMPILER_WORKS 1)
set(CMAKE_CXX_COMPILER_WORKS 1)

# Adjust the default behavior of the FIND_XXX() commands
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Default to Cortex-M4 if not specified
if(NOT DEFINED CORTEX_TYPE)
    set(CORTEX_TYPE "M4")
endif()

# Set CPU and FPU flags based on Cortex type
if(CORTEX_TYPE STREQUAL "M0")
    set(CPU_FLAGS "-mcpu=cortex-m0 -mthumb")
elseif(CORTEX_TYPE STREQUAL "M0PLUS")
    set(CPU_FLAGS "-mcpu=cortex-m0plus -mthumb")
elseif(CORTEX_TYPE STREQUAL "M3")
    set(CPU_FLAGS "-mcpu=cortex-m3 -mthumb")
elseif(CORTEX_TYPE STREQUAL "M4")
    set(CPU_FLAGS "-mcpu=cortex-m4 -mthumb -mfpu=fpv4-sp-d16 -mfloat-abi=hard")
elseif(CORTEX_TYPE STREQUAL "M7")
    set(CPU_FLAGS "-mcpu=cortex-m7 -mthumb -mfpu=fpv5-d16 -mfloat-abi=hard")
elseif(CORTEX_TYPE STREQUAL "M33")
    set(CPU_FLAGS "-mcpu=cortex-m33 -mthumb -mfpu=fpv5-sp-d16 -mfloat-abi=hard")
else()
    message(WARNING "Unknown CORTEX_TYPE: ${CORTEX_TYPE}, using default M4 flags")
    set(CPU_FLAGS "-mcpu=cortex-m4 -mthumb -mfpu=fpv4-sp-d16 -mfloat-abi=hard")
endif()

# Common flags for embedded targets
set(COMMON_FLAGS "${CPU_FLAGS} -fdata-sections -ffunction-sections")

set(CMAKE_C_FLAGS_INIT "${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS_INIT "${COMMON_FLAGS}")
set(CMAKE_ASM_FLAGS_INIT "${COMMON_FLAGS}")

# Linker flags for embedded targets
set(CMAKE_EXE_LINKER_FLAGS_INIT "${CPU_FLAGS} -Wl,--gc-sections -specs=nano.specs -specs=nosys.specs")
