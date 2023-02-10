#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "iter_e_ipdft_imp.h"

#define NUM_CHANNNELS 3

int main() {

    estimator_config pmu_config;
    
    //performance test
    clock_t start, end;
    unsigned int perf_iterations = 1000; 
    double avg_perf_time = 0;

    double AMP = 2;
    double PH = 0;
    double FREQ = 50;

    double ki = 0.1;
    double fi = 75;
    unsigned int n =2048 ;
    double fs = 25600;
    unsigned int n_bins = 11;    
    int P = 3;
    int Q = 22;
    double epsilon = 0.0033;
    double dt = 1/fs;
    double df = fs/n;

    int i, j;
    //allocating memory
    double amp[NUM_CHANNNELS];
    double ph[NUM_CHANNNELS];
    double freq[NUM_CHANNNELS]; 
    synchrophasor* estimated_phasors = (synchrophasor *)malloc(NUM_CHANNNELS * sizeof(synchrophasor));
    double** signal_windows = (double **)malloc(NUM_CHANNNELS * sizeof(double *));
    for (i = 0; i < NUM_CHANNNELS; i++){
        signal_windows[i] = (double *)malloc(n * sizeof(double));
    }
    //initializing windows
    for (j=0; j<NUM_CHANNNELS; j++){
        amp[j] = AMP;
        ph[j] = PH + j*2*M_PI/3;
        freq[j] = FREQ;
        for(i=0; i<n; i++){
            signal_windows[j][i] = (amp[j]*cos(2*M_PI*freq[j]*dt*i + ph[j]) + amp[j]*ki*cos(2*M_PI*fi*dt*i + ph[j]));
        }
    }

    printf("\n== M-Class Parameters ========================================================\n\n");
    for(j=0; j<NUM_CHANNNELS; j++){
       printf("[Channel:%d] Fundamental Component | Amp(V): %0.2lf | Ph(rad): %0.2lf | Freq(Hz): %0.2lf\n",j, amp[j], ph[j], freq[j]); 
    }
    printf("\nInterference | I-Mag(%%): %0.2lf | I-Freq(Hz): %0.2lf\n", ki*100, fi);
    printf("------------------------------------------------------------------------------\n");
    printf("Window | SamplingFreq(kS/s): %0.3lf | NCycles: %1.0f | FreqResolution: %0.2lf\n", (float)fs/1000, (n*50/fs), df);
    printf("Iterations | P: %d | Q: %d \n", P, Q);
    printf("\n===============================================================================\n");

    pmu_config.win_len = n;
    pmu_config.fs = fs;
    pmu_config.n_bins = n_bins;
    pmu_config.P = P;
    pmu_config.Q = Q;
    pmu_config.interf_trig = epsilon;
    pmu_config.n_chanls = NUM_CHANNNELS;

    pmu_init(&pmu_config);
    for (i= 0; i<perf_iterations; i++){
        start = clock();
        pmu_estimate(signal_windows, estimated_phasors);
        end = clock();
        avg_perf_time += (double)(end - start) / CLOCKS_PER_SEC;
    }
    avg_perf_time = avg_perf_time/perf_iterations;

    printf("\n---- [Results]: -------------------------------------------------------------------------------------\n\n");
    printf("| Number of Itrations: %d \n", perf_iterations);
    printf("| Total Estimation Time (seconds): %.10lf \n", avg_perf_time);
    for(j=0; j<NUM_CHANNNELS; j++){
        printf("| CHANNEL: %d |\tFREQ: %.10lf (Hertz) | AMP: %.5lf (Volt) | PH: %.10lf (deg)\n",j, estimated_phasors[j].freq, estimated_phasors[j].amp, estimated_phasors[j].ph*(180/M_PI));
    }
    printf("\n-----------------------------------------------------------------------------------------------------\n");

    pmu_deinit();

    for (i = 0; i < NUM_CHANNNELS; i++){
        free(signal_windows[i]);
    }
    free(signal_windows);
    free(estimated_phasors);

    return 0;
}
