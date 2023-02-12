# __m-class-pmu__
An ANSI C implementation of the Phasor Measurment Unit Estimator (PMU Estimator) based on the Iterative Interpolated DFT Synchrophasor Estimation Algorithm.
# Version 1.4.1
Updates:

- Now the library supprots CMake Building!
- The inlined functions are now defined as preprocessor macros.


## __Building the library__
To build the library, first make sure that you have the following build tools are installed:

- CMake (VERSION 3.0 or higher)
- GCC

You can build the library either using CMake or by just running the __build.sh__ provided in the root directory of the project.
### __Building using build.sh__

run the following command from the root directory of the project:

    (Linux: sudo chmod +x build.sh)
    ./build.sh

This will build both the static and shared libraries. The libraries will be placed in the __/build__ directory along with the __pmu_estimator.h__ header and the __config.ini__ estimator configuration file.

### __Building using CMake__

run the following commands from the root directory of the project:

    mkdir build
    cd build
    cmake ..

then to build as a static library run:

    make PmuEstimatorStatic

or to build as a shared library run:

    make PmuEstimatorShared

the libraries will be placed in the __/build__.

## __Pmu Estimator Configuration__

In order to configure the program, the __config.ini__ which can be found in __config/config.ini__ file must be edited. otherwise pass config structure to the __pmu_init()__ function.

To specify whether the program should use the config structure passed as argument or the one specified in the __config.ini__ file do as follows:

for the __config.ini__ file:

    char filename[] = "config/config.ini";
    pmu_init(&filename, CONFIG_FROM_INI);

for the config structure:
    
    pmu_init(&config, CONFIG_FROM_STRUCT);

## __Enabling Debug Logs__
In order to compile the library with debug logs enabled, the __DEBUG__ directive must be defined. To do so, you can use the -D flag when runnin the ./build.sh command as follows:

    ./build.sh -D

or add the __-DDEBUG=ON__ with the __cmake__:

    cmake -DDEBUG=ON ..
## __Running the "main" _program (tested on Windows)__
to compile and run the code simply run:

    gcc -I .\libs\iniparser -I .\src .\libs\iniparser\dictionary.c .\libs\iniparser\iniparser.c main.c .\src\pmu_estimator.c -o main.exe
    ./main.exe

## __Computational Time Evaluation__
To keep track of the computational time per call of the __pmu_estimate()__ function for a certain number of calls, a timer is set.
To change the number of calls to perform in order to compute the average time per call set the value of the following directive:

    #define PERF_ITERATIONS 1000


