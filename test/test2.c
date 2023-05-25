#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pmu_estimator.h"

//#define DEBUG_LOG

#define NUM_CHANNNELS 1
#define PERF_ITERATIONS 1000
#define MAX_N_CYCLES 16
#define MAX_Q 50
#define MAX_FS 51200
#define MIN_N_BINS 11

#define wrap_angle(rad) (((double)(rad) - 2 * M_PI * rint((double)(rad) / (2 * M_PI)))) // wraps angle (rad) in range [-pi; pi]

int main() {

    // Creating a file to store the results =================================================

    FILE *file;
    char *filename = "pmu_perf.csv";
    file = fopen(filename, "w");
    if (file == NULL) {
        printf("Failed to create the file.\n");
        return 1;
    }
    fprintf(file, "Sample Rate (S/s), Number of Cycles, Q Iterations (OOBI), Number of DFT bins, avg pmu_estimate() wall time (ms), Estimated Freq, Estimated Amp, Estimated Ph, Estimated Rocof, Expected Freq, Expected Amp, Expected Ph, Expected Rocof \n");

    estimator_config pmu_config;
    
    //performance test
    struct timespec start_ns, end_ns;

    double AMP = 2;
    double PH = 0;
    double FREQ = 51;

    unsigned int f0 = 50;
    unsigned int frame_rate = 50 ;
    double ki = 0.1;
    double fi = 75;

    _Bool iter_eipdft = 1;  
    int P = 3;
    
    double epsilon = 0.0033;
    double th_coeff[3]= {3, 25, 0.035};
    double lpf_coeff[3]= {0.5913, 0.2043, 0.2043};
    
    int total_iterations = MAX_FS/25600 * ((MAX_N_CYCLES)/2) * (MAX_Q/10) * 2;
    int iter_counter = 0;

    for (unsigned int fs=25600; fs <= MAX_FS; fs=fs*2){

        for (int n_cycles=2; n_cycles <= MAX_N_CYCLES; n_cycles=n_cycles*2){

            double total_perf_time_ns = 0;
            double avg_perf_time_ns = 0;

            for (int Q = 3 ; Q <= MAX_Q; Q = Q+10){

                unsigned int n_bins = (n_cycles + 2 < MIN_N_BINS) ? MIN_N_BINS : n_cycles + 2;

                double n = (double)n_cycles*fs/f0;
                double dt = 1/(double)fs;
                double df = (double)fs/n;

                int i, j, chanl;
                //allocating memory
                double amp[NUM_CHANNNELS];
                double ph[NUM_CHANNNELS];
                double freq[NUM_CHANNNELS]; 
                pmu_frame estimated_frame[NUM_CHANNNELS];
                float_p signal_windows[NUM_CHANNNELS][(int)n];

                //initializing signal windows
                for (chanl=0; chanl<NUM_CHANNNELS; chanl++){
                    
                    amp[chanl] = AMP;
                    ph[chanl] = PH + chanl*2*M_PI/3;
                    freq[chanl] = FREQ;
                    for(j=0; j<n; j++){
                        signal_windows[chanl][j] = (amp[chanl]*cos(2*M_PI*freq[chanl]*dt*j + ph[chanl]) + amp[chanl]*ki*cos(2*M_PI*fi*dt*j + ph[chanl]));
                    }
                }

#ifdef DEBUG_LOG
                printf("\n== PMU-Estimator ============================================================\n\n");
                for(j=0; j<NUM_CHANNNELS; j++){
                printf("[Channel:%d] Fundamental Component | Amp(V): %0.2lf | Ph(rad): %0.2lf | Freq(Hz): %0.2lf\n",j, amp[j], ph[j], freq[j]); 
                }
                printf("\nInterference | I-Mag(%%): %0.2lf | I-Freq(Hz): %0.2lf\n", ki*100, fi);
                printf("------------------------------------------------------------------------------\n");
                printf("Window | SamplingFreq(kS/s): %0.3lf | NCycles: %u | FreqResolution: %0.2lf\n", (float)fs/1000, n_cycles, df);
                printf("Iterations | P: %d | Q: %d \n", P, Q);
                printf("\n===============================================================================\n");           
#endif
                // Configuring PMU Estimator =========================================================

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

                // Computation Iterations =============================================================

                double mid_fracsec = 0;
                for (i= 0; i<PERF_ITERATIONS; i++){
                    timespec_get (&start_ns, TIME_UTC);

                    pmu_estimate((float_p *)signal_windows, mid_fracsec ,estimated_frame);
                    
                    timespec_get (&end_ns, TIME_UTC);
                    total_perf_time_ns += (end_ns.tv_sec - start_ns.tv_sec) + (end_ns.tv_nsec - start_ns.tv_nsec) / 1000000000.0;
                }
                avg_perf_time_ns = total_perf_time_ns/PERF_ITERATIONS;


                // Results ============================================================================
#ifdef DEBUG_LOG
                printf("\n---- [Results]: -------------------------------------------------------------------------------------\n\n");
                printf("| Number of Itrations: %d \n", PERF_ITERATIONS);

                printf("| HP: Total Estimation Time (ms): %.10lf \n", 1000*total_perf_time_ns);
                printf("| HP: Avg Estimation Time per Call (ms): %.10lf \n\n", 1000*avg_perf_time_ns);

                for(j=0; j<NUM_CHANNNELS; j++){

                    estimated_frame[j].synchrophasor.ph = wrap_angle(estimated_frame[j].synchrophasor.ph);
                    printf("| CHANNEL: %d |\tFREQ: %.10lf (Hertz) | AMP: %.5lf (Volt) | PH: %.10lf (deg) | ROCOF: %.10lf (Hz/s)\n",j, estimated_frame[j].synchrophasor.freq, estimated_frame[j].synchrophasor.amp, estimated_frame[j].synchrophasor.ph*(180/M_PI), estimated_frame[j].rocof);
                }
                printf("\n-----------------------------------------------------------------------------------------------------\n");
#endif
                // Writing the results to the file =====================================================
                fprintf(file, "%d, %d, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", fs, n_cycles, Q, n_bins, 1000*avg_perf_time_ns ,estimated_frame[0].synchrophasor.freq, estimated_frame[0].synchrophasor.amp, estimated_frame[0].synchrophasor.ph*(180/M_PI), estimated_frame[0].rocof, freq[0], amp[0], ph[0]*(180/M_PI), 2*M_PI*f0);

                iter_counter ++;
                printf("Progress: %.4lf %%\n", 100*iter_counter/(double)total_iterations);
                
                // Deinitializing PMU Estimator =========================================================
                pmu_deinit();

            } 
        }
    }

    fclose(file);
    
    return 0;
}
