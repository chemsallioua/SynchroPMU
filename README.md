# __SynchroPMU__
A __C implementation__ of the Phasor Measurment Unit Estimator (PMU Estimator) based on the Iterative Interpolated DFT Synchrophasor Estimation Algorithm.
# Version 1.7.0
Updates (with respect to version 1.6.4):

- Now it is possible to create independent instances of the pmu estimator which have their own configuration and state. This works both in the native C library and in the Python API.


previous updates (with respect to version 1.3.0):
- Performance evaluation by means of __test2.c__ added to the __/test__ directory.
- Python API now available (ver:1.6.3)! check in _/python_
- updated the library license (BSD 3-clause) 
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
- now it's easily possible to stub the arithmetic functions with user implementations.
- added a test script and a Makefile to test the library with __gprof__.
- optimized the FFT function to work faster when input window is power of 2, and also added the -03 flag to GCC for optimized compiling.
- the library supports logging with different levels: 0 (no logging), 1 (only errors), 2 (errors and info logs), 3 (errors, info and debug logs).
- further optimized the __pmu_estimate()__ function by optimizing wf() hanning FT function.

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

__NOTE__: if you set the number of channels different from the dimention of the input signal array passed to the __pmu_estimate()__ function, you experience  ___Segmentation Error___.

## __Enabling Logging__
To compile the library with logs enabled, the __LOGGING_LEVEL__ directive must be defined and set to a logging level. To do so, you can use the -D option foollowed by the logging level wanted, when running the ./build.sh. 

This are the available logging levels:
    
- __LOGGING_LEVEL = 0__ : No logging
- __LOGGING_LEVEL = 1__ : Only errors
- __LOGGING_LEVEL = 2__ : Errors and Info logs
- __LOGGING_LEVEL = 3__ : Errors, Info and Debug logs (which include also each step taken by the algorithm and the various algorithm variable values)

For example, to build with logging level 3, run:

    ./build.sh -D 3

or if you are building directly with __cmake__ add the __-DLOGGING_LEVEL=3__ option:

    cmake -DLOGGING_LEVEL=3 ..

## __Stub Arithmetic functions__
__func_stub.c__ is a header file that defines macros to stub out several data types and functions.
 In the file you can use the macros to either implement your own version of a function or to replace it with a version from an external library.
 
 __For example:__ if you want to implement your own version of the __pmue_fft_r__ function, you can define a new function with the same name and signature, and then use the __pmue_fft_r__ macro to replace the original function with your new implementation:

    #define pmue_fft_r(in_ptr, out_ptr, out_len, n_bins) my_fft((in_ptr), (out_ptr), (out_len), (n_bins))

    int_p my_fft(float_p* in_ptr, float_p complex_p* out_ptr , uint_p out_len, uint_p n_bins) {
        //Your implementation of fft here
    }

Note that if you want to use a function from an external library, you will need to link the library to your code during the compilation process.
## __Pmu Estimator Configuration__

In order to configure the program, the __config.ini__ which can be found in __config/config.ini__ file must be edited. otherwise pass config structure to the __pmu_init()__ function.

To specify whether the program should use the config structure passed as argument or the one specified in the __config.ini__ file do as follows:

for the __config.ini__ file:

    char filename[] = "config/config.ini";
    pmu_init(&filename, CONFIG_FROM_INI);

for the config structure:
    
    pmu_init(&config, CONFIG_FROM_STRUCT);

## __Building and Running the examples__

To build an example found in the __examples__ directory, you can build it by running the following command from the __/examples__ directory of the project:

    make <name of the example>

or to build all examples simply run:

    make

the executable will be placed in the __/examples/build__ directory.
then run the executable:
    
    ./build/<name of the example>

## __Testing, Profiling and Performance Evaluation__

### __Profiling with gprof__
In the folder __test__ you can find a test script that can be used to test and profile library with __gprof__. To run the test script, first change the configuration of the pmu estimator in the __test.c__ file and also the number of test iterations for the __pmu_estimate()__ function by setting the __PERF_ITERATIONS__ directive constant, then run from the __/test__ directory:

    make test

this will build the test executable, now run:
    
    make profile

this will run the test executable and generate a  __profile.txt__ which contain the profiling results. To clean the test directory run:

    make clean

### __Performance Evaluation__
In order to evaluate the computation time in other words the wall time for the excecution of the pmu_estimate() function, you can customize the file __test2.c__ in the __/test__ directory. In the file you can set the number of iterations for the __pmu_estimate()__ function by setting the __PERF_ITERATIONS__ directive constant and other parameters as well. Then run the following command from the __/test__ directory:

    make test2

this will build the test executable, now run:
        
    make benchmark

this will run the test executable and generate a  __pmu_perf.csv__ which contain the results. To clean the test directory run:

    make clean


# __Python API (v1.2.0)__

The library also comes with a python API that can be used to call the __pmu_estimate()__ function from python. To use the API, first make sure that you have the following build tools are installed:

- CMake (VERSION 3.0 or higher)
- A C compiler, such as GCC
- Make
- Python3 (with ctypes)

Make sure that you have already built the library using the __build.sh__ script or using __cmake__ and also installed it, specifically the dynamic library. (See __Building the library__ section for more details.)

Now install the library from source by running the following command from the root directory of the project:

    pip install ./python

This will install the library as a python package. Now you can import the package and use it in your python code.

## __Running the exaample__

After installing the python package, you can take a loot and run the example found in the __/python/examples__ directory by running the following command from the __/python/examples__ directory:

    python ./python/examples/example_usage.py
