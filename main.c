#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "iter_e_ipdft_imp.h"

int main() {

    estimator_config pmu_config;

    double amp=2;
    double ph=1;
    double freq= 50;

    double ki = 0.1;
    double fi = 25;
    double n =2048 ;
    double fs = 25600;
    unsigned int n_bins = 11;    
    int P = 3;
    int Q = 3;
    double epsilon = 0.0033;

    double signal_window[(int)n];

    double dt = 1/fs;
    double df = fs/n;

    int i, j;
    for(i=0; i<n; i++){
        signal_window[i] = (amp*cos(2*M_PI*freq*dt*i + ph) + amp*ki*cos(2*M_PI*fi*dt*i + ph));
    }

    pmu_config.win_len = n;
    pmu_config.fs = fs;
    pmu_config.n_bins = n_bins;
    pmu_config.P = P;
    pmu_config.Q = Q;
    pmu_config.interf_trig = epsilon;

    printf("\n== M-Class Parameters ========================================================\n");
    printf("Signal Fundamental Component | Amp(V): %0.2lf | Ph(rad): %0.2lf | Freq(Hz): %0.2lf\n", amp, ph, freq);
    printf("Interference | I-Mag(%%): %0.2lf | I-Freq(Hz): %0.2lf\n", ki*100, fi);
    printf("------------------------------------------------------------------------------\n");
    printf("Window | SamplingFreq(kS/s): %0.3lf | NCycles: %1.0f | FreqResolution: %0.2lf\n", (float)fs/1000, (n*50/fs), df);
    printf("Iterations | P: %d | Q: %d \n", P, Q);
    printf("===============================================================================\n");

    pmu_init(&pmu_config);

    synchrophasor estimated_phasor;
    pmu_estimate(signal_window, &estimated_phasor);

    printf("\n---- [Results] ----------------------------------------------------------\n");
    printf("|\tFREQ: %.10lf | AMP: %.10lf | PH: %.10lf \t|\n", estimated_phasor.freq, estimated_phasor.amp, estimated_phasor.ph);
    printf("-------------------------------------------------------------------------\n");

    pmu_deinit();

    return 0;
}
