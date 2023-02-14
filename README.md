# PMU Estimator
An ANSI C implementation of the Phasor Measurment Unit Estimator (PMU Estimator) based on the Iterative Interpolated DFT Synchrophasor Estimation Algorithm.
# Version 1.4.8
Updates (with respect to version 1.3.0):

- Now the library supprots CMake Building!
- The inlined functions are now defined as preprocessor macros.
- the __pmu_estimate()__ has an additional input argument, __mid_fracsec__ which is the fraction of second relative to PPS of the mid point of the window, it is used t make the estimated phase correction.
- added a flag in the pmu configuration structure to specify whether the enhanced iterative interpolated dft should be applied or not.
- added two new config files aimed at specifying specific configurations for M-Class pmu and for P-Class pmu.
- fixed a bug in the __pmu_estimate()__ function that caused the program to raise an error when the input signal window is passed as a static array.
- added two examples that shows how to use the library with __CONFIG_FROM_STRUCT__ and also with __CONFIG_FROM_INI__.
- fixed bug of __wrap_angle()__ low result accuracy.
- fixed bug in CmakeLists.txt that caused the library to raise an error when building with __NUM_CHANLS__ not set.
- added library installation with cmake.

## __Building the library__
To build the library, first make sure that you have the following build tools are installed:

- CMake (VERSION 3.0 or higher)
- A C compiler, such as GCC
- Make

You can build the library either using CMake or by just running the __build.sh__ provided in the root directory of the project.
### __Building using build.sh__

run the following command from the root directory of the project:

    (Linux: sudo chmod +x build.sh)
    ./build.sh

This will build both the static and shared libraries. The libraries will be placed in the __/build__ directory along with the __pmu_estimator.h__ header and the __config.ini__ estimator configuration file. This will also install the libraries as it executes the __cmake --install .__ command.

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

### __Installing the library__
use the following command to install the library from the __/build__ directory:

    cmake --install .
## __Setting Number of Channels__
To set number of channels on which the pmu estimator will process with frames, __NUM_CHANLS__ directive must be defined. The default value is __NUM_CHANLS = 1__. To set the value, you can use the -N option when running the ./build.sh command.for example, to set the number of channels to 4:

    ./build.sh -N 4

or add the __-DNUM_CHANLS=4__ with the __cmake__:

    cmake -DNUM_CHANLS=4 ..

## __Enabling Debug Logs__
To compile the library with debug logs enabled, the __DEBUG__ directive must be defined. To do so, you can use the -D flag when running the ./build.sh command as follows:

    ./build.sh -D

or add the __-DDEBUG=ON__ with the __cmake__:

    cmake -DDEBUG=ON ..

## __Pmu Estimator Configuration__

In order to configure the program, the __config.ini__ which can be found in __config/config.ini__ file must be edited. otherwise pass config structure to the __pmu_init()__ function.

To specify whether the program should use the config structure passed as argument or the one specified in the __config.ini__ file do as follows:

for the __config.ini__ file:

    char filename[] = "config/config.ini";
    pmu_init(&filename, CONFIG_FROM_INI);

for the config structure:
    
    pmu_init(&config, CONFIG_FROM_STRUCT);

## __Building and Running the examples__

To run an example found in the __examples__ directory, you can build it by running the following command from the root directory of the project:

    gcc -I ./libs/iniparser -I ./src ./libs/iniparser/dictionary.c ./libs/iniparser/iniparser.c ./examples/<name of the example c file> ./src/pmu_estimator.c -o <name of the example>

then run the executable:
    
        ./<name of the example>
