# __m-class-pmu__
A C implementtation of the M-Class Iterative Interpolated DFT Synchrophasor Estimattion Algorithm
# Version 1.0.0
Updates:
- Added debug logging feature
- Adapted the algorithm to match the Iterative Enhanced Interpolated DFT

## __Activating Debug Logs__
change the value of the _#define DEBUG 0_ in the __Mclass.c__ file

to deactivate:

    #define DEBUG 0
or to activate

    #define DEBUG 1

## __Running the program (tested on Windows)__
to compile and run the code simply run:

    gcc .\Mclass.c -o Mclass | ./Mclass.exe