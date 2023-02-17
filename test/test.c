#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pmu_estimator.h"

#define NUM_CHANNNELS 1
#define PERF_ITERATIONS 10000

int main() {

    estimator_config pmu_config;
    
    //performance test
    clock_t start, end;  
    double avg_perf_time = 0;

    double AMP = 2;
    double PH = 0;
    double FREQ = 51;

    unsigned int f0 = 50;
    unsigned int frame_rate = 50 ;
    double ki = 0;
    double fi = 75;
    unsigned int n_cycles = 4;
    unsigned int fs = 25600;
    unsigned int n_bins = 6;  
    _Bool iter_eipdft = 1;  
    int P = 3;
    int Q = 22;
    double epsilon = 0.0033;
    double th_coeff[3]= {3, 25, 0.035};
    double lpf_coeff[3]= {0.5913, 0.2043, 0.2043};


    double n = (double)n_cycles*fs/f0;
    double dt = 1/(double)fs;
    double df = (double)fs/n;

    int i, j, chanl;
    //allocating memory
    double amp[NUM_CHANNNELS];
    double ph[NUM_CHANNNELS];
    double freq[NUM_CHANNNELS]; 
    pmu_frame estimated_frame[NUM_CHANNNELS];
    double signal_windows[NUM_CHANNNELS][(int)n];

    //printf("channels: %d, window_size: %u\n", NUM_CHANNNELS, (unsigned int)n);

    //initializing windows
    for (chanl=0; chanl<NUM_CHANNNELS; chanl++){
        amp[chanl] = AMP;
        ph[chanl] = PH + chanl*2*M_PI/3;
        freq[chanl] = FREQ;
        for(j=0; j<n; j++){
            signal_windows[chanl][j] = (amp[chanl]*cos(2*M_PI*freq[chanl]*dt*j + ph[chanl]) + amp[chanl]*ki*cos(2*M_PI*fi*dt*j + ph[chanl]));
        }
    }

    pmu_config.n_cycles = n_cycles;
    pmu_config.f0 = f0;
    pmu_config.frame_rate = frame_rate;
    pmu_config.fs = fs;
    pmu_config.n_bins = n_bins;
    pmu_config.P = P;
    pmu_config.Q = Q;
    pmu_config.interf_trig = epsilon;
    pmu_config.iter_eipdft = iter_eipdft;

    pmu_config.rocof_thresh[0] = th_coeff[0];
    pmu_config.rocof_thresh[1] = th_coeff[1];
    pmu_config.rocof_thresh[2] = th_coeff[2];

    pmu_config.rocof_low_pass_coeffs[0] = lpf_coeff[0];
    pmu_config.rocof_low_pass_coeffs[1] = lpf_coeff[1];
    pmu_config.rocof_low_pass_coeffs[2] = lpf_coeff[2];

    pmu_init(&pmu_config, CONFIG_FROM_STRUCT);

    double mid_fracsec = 0;
    for (i= 0; i<PERF_ITERATIONS; i++){
        start = clock();

        pmu_estimate((double *)signal_windows, mid_fracsec ,estimated_frame);

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
