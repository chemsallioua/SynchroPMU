# __m-class-pmu__
A C implementtation of the M-Class Iterative Interpolated DFT Synchrophasor Estimattion Algorithm
# Version 1.3.0
Updates:

- Possibility of configuring the pmu estimator via both INI file or programmatically using the config struct
- Validity check of the config values
- reorganized the program architecture
- pmu estimator uses /iniparser C library to parse the config INI file.

## __Activating Debug Logs__
change the value of the _#define DEBUG 0_ in the __iter_e_ipdft_imp.h__ file

to deactivate:

    #define DEBUG 0
or to activate

    #define DEBUG 1

## __Running the program (tested on Windows)__
to compile and run the code simply run:

    gcc -I .\libs\iniparser .\libs\iniparser\dictionary.c .\libs\iniparser\iniparser.c main.c pmu_estimator.c -o main.exe
    ./main.exe

## __Computational Time Evaluation__
To keep track of the computational time per call of the __pmu_estimate()__ function for a certain number of calls, a timer is set.
To change the number of calls to perform in order to compute the average time per call set the value of the following directive:

    #define PERF_ITERATIONS 1000

## __Configuration__

In order to configure the program, the __config.ini__ which can be found in __config/config.ini__ file must be edited. otherwise pass config structure to the __pmu_init()__ function.

To specify whether the program should use the config structure passed as argument or the one specified in the __config.ini__ file do as follows:

for the __config.ini__ file:

    char filename[] = "config/config.ini";
    pmu_init(&filename, CONFIG_FROM_INI);

for the config structure:
    
    pmu_init(&config, CONFIG_FROM_STRUCT);
