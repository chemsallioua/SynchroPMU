# __m-class-pmu__
A C implementtation of the M-Class Iterative Interpolated DFT Synchrophasor Estimattion Algorithm
# Version 1.1.1
Updates:
- feature: computation on multiple configurable channels
## __Activating Debug Logs__
change the value of the _#define DEBUG 0_ in the __iter_e_ipdft_imp.h__ file

to deactivate:

    #define DEBUG 0
or to activate

    #define DEBUG 1

## __Running the program (tested on Windows)__
to compile and run the code simply run:

    gcc main.c iter_e_ipdft_imp.c -o main.exe
    ./main.exe
