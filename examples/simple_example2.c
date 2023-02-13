#include <stdio.h>
#include <math.h>

#include "pmu_estimator.h"

int main() {

    // Input Signal parameters
    double AMP = 2;
    double PH = 1;
    double FREQ = 50;
    unsigned int sample_rate = 25600;
    unsigned int window_size = 2048;
    double dt = 1/(double)sample_rate;

    double input_signal_window[window_size];

    // Initializing input signal window
    unsigned int i;
    for(i=0; i<window_size; i++){
        input_signal_window[i] = AMP*cos(2*M_PI*FREQ*dt*i + PH);
    }

    // Initializing pmu estimator from config struct
    estimator_config pmu_config;

    pmu_config.n_cycles = 4;
    pmu_config.f0 = 50;
    pmu_config.frame_rate = 50;
    pmu_config.fs = 25600;
    pmu_config.n_bins = 11;
    pmu_config.P = 3;
    pmu_config.Q = 3;
    pmu_config.interf_trig = 0.0033;
    pmu_config.iter_eipdft = 1;

    pmu_config.rocof_thresh[0] = 3;
    pmu_config.rocof_thresh[1] = 25;
    pmu_config.rocof_thresh[2] = 0.035;

    pmu_config.rocof_low_pass_coeffs[0] = 0.5913;
    pmu_config.rocof_low_pass_coeffs[1] = 0.2043;
    pmu_config.rocof_low_pass_coeffs[2] = 0.2043;

    pmu_init(&pmu_config, CONFIG_FROM_STRUCT);

    pmu_frame estimated_frame;
    double mid_window_fracsec = 0;

    // Estimating frame
    pmu_estimate((double (*)[1])input_signal_window, mid_window_fracsec , &estimated_frame);
    
    // Printing estimated frame
    pmu_dump_frame(&estimated_frame, stdout);

    // Deinitializing pmu estimator
    pmu_deinit();

    return 0;
}