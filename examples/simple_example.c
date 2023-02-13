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

    // Initializing pmu estimator from config_example.ini
    char file_name[] = "examples/config_example.ini";
    pmu_init(&file_name, CONFIG_FROM_INI);

    pmu_frame estimated_frame;
    double mid_window_fracsec = 0;

    // Estimating frame
    pmu_estimate((double *)input_signal_window, mid_window_fracsec , &estimated_frame);
    
    // Printing estimated frame
    pmu_dump_frame(&estimated_frame, stdout);

    // Deinitializing pmu estimator
    pmu_deinit();

    return 0;
}