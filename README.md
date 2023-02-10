# __m-class-pmu__
A C implementtation of the M-Class Iterative Interpolated DFT Synchrophasor Estimattion Algorithm
# Version 1.2.0
Updates:
- new feature: now the pmu_estimate() function computes also the ROCOF
## __Activating Debug Logs__
change the value of the _#define DEBUG 0_ in the __iter_e_ipdft_imp.h__ file

to deactivate:

    #define DEBUG 0
or to activate

    #define DEBUG 1

## __Running the program (tested on Windows)__
to compile and run the code simply run:

    gcc main.c iter_e_ipdft_imp.c -o main.exe | ./main.exe

## __Computational Time Evaluation__
To keep track of the computational time per call of the __pmu_estimate()__ function for a certain number of calls, a timer is set.
To change the number of calls to perform in order to compute the average time per call set the value of the following directive:

    #define PERF_ITERATIONS 1000