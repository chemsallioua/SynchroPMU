#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pmu_estimator.h"

#define NUM_CHANNNELS 6
#define PERF_ITERATIONS 1000

int main() {

    estimator_config pmu_config;
    
    //performance test
    clock_t start, end;  
    double avg_perf_time = 0;

    double AMP = 2;
    double PH = 0;
    double FREQ = 50;

    unsigned int f0 = 50;
    unsigned int frame_rate = 50 ;
    double ki = 0.1;
    double fi = 75;
    unsigned int n_cycles = 4;
    unsigned int fs = 25600;
    unsigned int n_bins = 11;  
    _Bool iter_eipdft = 1;  
    int P = 3;
    int Q = 22;
    double epsilon = 0.0033;
    double th_coeff[3]= {3, 25, 0.035};
    double lpf_coeff[3]= {0.5913, 0.2043, 0.2043};


    double n = (double)n_cycles*fs/f0;
    double dt = 1/(double)fs;
    double df = (double)fs/n;

    int i, j;
    //allocating memory
    double amp[NUM_CHANNNELS];
    double ph[NUM_CHANNNELS];
    double freq[NUM_CHANNNELS]; 
    pmu_frame estimated_frame[NUM_CHANNNELS];
    double signal_windows[NUM_CHANNNELS][(int)n];

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
    printf("Window | SamplingFreq(kS/s): %0.3lf | NCycles: %u | FreqResolution: %0.2lf\n", (float)fs/1000, n_cycles, df);
    printf("Iterations | P: %d | Q: %d \n", P, Q);
    printf("\n===============================================================================\n");

    pmu_config.n_cycles = n_cycles;
    pmu_config.f0 = f0;
    pmu_config.frame_rate = frame_rate;
    pmu_config.fs = fs;
    pmu_config.n_bins = n_bins;
    pmu_config.P = P;
    pmu_config.Q = Q;
    pmu_config.interf_trig = epsilon;
    pmu_config.n_chanls = NUM_CHANNNELS;
    pmu_config.iter_eipdft = iter_eipdft;

    pmu_config.rocof_thresh[0] = th_coeff[0];
    pmu_config.rocof_thresh[1] = th_coeff[1];
    pmu_config.rocof_thresh[2] = th_coeff[2];

    pmu_config.rocof_low_pass_coeffs[0] = lpf_coeff[0];
    pmu_config.rocof_low_pass_coeffs[1] = lpf_coeff[1];
    pmu_config.rocof_low_pass_coeffs[2] = lpf_coeff[2];

    //pmu_init(&pmu_config, CONFIG_FROM_STRUCT);
    char file_name[] = "config/config.ini";
    pmu_init(&file_name, CONFIG_FROM_INI);
    double mid_fracsec = 0;
    for (i= 0; i<PERF_ITERATIONS; i++){
        start = clock();

        if(pmu_estimate(signal_windows, mid_fracsec ,estimated_frame)){
            printf("Error: Estimation Failed!\n");
            return -1;
        }
        end = clock();
        avg_perf_time += (double)(end - start) / CLOCKS_PER_SEC;
    }
    avg_perf_time = avg_perf_time/PERF_ITERATIONS;

    printf("\n---- [Results]: -------------------------------------------------------------------------------------\n\n");
    printf("| Number of Itrations: %d \n", PERF_ITERATIONS);
    printf("| Total Estimation Time (seconds): %.10lf \n", avg_perf_time);
    for(j=0; j<NUM_CHANNNELS; j++){
        printf("| CHANNEL: %d |\tFREQ: %.10lf (Hertz) | AMP: %.5lf (Volt) | PH: %.10lf (deg) | ROCOF: %.10lf (Hz/s)\n",j, estimated_frame[j].synchrophasor.freq, estimated_frame[j].synchrophasor.amp, estimated_frame[j].synchrophasor.ph*(180/M_PI), estimated_frame[j].rocof);
    }
    printf("\n-----------------------------------------------------------------------------------------------------\n");

    pmu_deinit();

    return 0;
}
