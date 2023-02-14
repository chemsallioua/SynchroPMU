cmake_minimum_required(VERSION 3.0)

project(PmuEstimator)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-g -ansi -fPIC -Wall -Wextra -Werror -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function")

option(DEBUG "Build with debug logs" OFF)
set(NUM_CHANLS 1 CACHE STRING "Number of channels")

message(STATUS "NUM_CHANLS=${NUM_CHANLS}")
add_definitions(-DNUM_CHANLS=${NUM_CHANLS})

if (DEBUG)
    message(STATUS "DEBUG=1")
    add_definitions(-DDEBUG=1)
else ()
    message(STATUS "DEBUG=0")
    add_definitions(-DDEBUG=0)
endif ()

# Add source files
set(PMUESTIMATOR_SRCS
        src/pmu_estimator.c
        libs/iniparser/dictionary.c
        libs/iniparser/iniparser.c
        )

# Add header files
set(PMUESTIMATOR_HEADERS
        src/pmu_estimator.h
        libs/iniparser/dictionary.h
        libs/iniparser/iniparser.h
        )

include_directories(libs/iniparser)

# Create the static library
add_library(PmuEstimatorStatic STATIC ${PMUESTIMATOR_SRCS} ${PMUESTIMATOR_HEADERS})
set_target_properties(PmuEstimatorStatic PROPERTIES OUTPUT_NAME "pmu_estimator")

# Create the shared library
add_library(PmuEstimatorShared SHARED ${PMUESTIMATOR_SRCS} ${PMUESTIMATOR_HEADERS})
set_target_properties(PmuEstimatorShared PROPERTIES OUTPUT_NAME "pmu_estimator")
set_target_properties(PmuEstimatorShared PROPERTIES SOVERSION 1)

# Add an alias for the shared library so that other targets can depend on the shared library
add_library(PmuEstimator::Shared ALIAS PmuEstimatorShared)

# Add an include directory
target_include_directories(PmuEstimatorStatic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(PmuEstimatorShared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Install the libraries and header files
install(TARGETS PmuEstimatorStatic PmuEstimatorShared
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

install(FILES ${PMUESTIMATOR_HEADERS} DESTINATION include)

# Copy the config.ini file to the build directory
configure_file(config/config.ini config.ini COPYONLY)
configure_file(config/m_class_config.ini m_class_config.ini COPYONLY)
configure_file(config/p_class_config.ini p_class_config.ini COPYONLY)
configure_file(src/pmu_estimator.h ${CMAKE_BINARY_DIR}/pmu_estimator.h COPYONLY)
