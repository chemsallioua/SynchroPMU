#include <stdio.h>
#include <math.h>

#include "pmu_estimator.h"

int main() {

    // Input Signal parameters
    float_p AMP = 2;
    float_p PH = 1;
    float_p FREQ = 50;
    uint_p sample_rate = 25600;
    uint_p window_size = 2048;
    float_p dt = 1/(float_p)sample_rate;

    float_p input_signal_window[window_size];

    // Initializing input signal window
    uint_p i;
    for(i=0; i<window_size; i++){
        input_signal_window[i] = AMP*cos(2*M_PI*FREQ*dt*i + PH);
    }

    // Initializing pmu estimator from config_example.ini
    char file_name[] = "examples/config_example.ini";
    pmu_init(&file_name, CONFIG_FROM_INI);

    pmu_frame estimated_frame;
    float_p mid_window_fracsec = 0;

    // Estimating frame
    pmu_estimate((float_p *)input_signal_window, mid_window_fracsec , &estimated_frame);
    
    // Printing estimated frame
    pmu_dump_frame(&estimated_frame, stdout);

    // Deinitializing pmu estimator
    pmu_deinit();

    return 0;
}