#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "pmu_estimator.h"

int main() {

    double AMP = 2;
    double PH = 0;
    double FREQ = 50;
    unsigned int fs = 25600;
    unsigned int n = 2048;
    double dt = 1/(double)fs;

    double* input_signal_window = (double *)malloc(n * sizeof(double));

    // Initializing input signal window
    unsigned int i;
    for(i=0; i<n; i++){
        input_signal_window[j][i] = AMP*cos(2*M_PI*FREQ*dt*i + PH);
    }

    // Initializing pmu estimator
    char file_name[] = "config/config.ini";
    pmu_init(&file_name, CONFIG_FROM_INI);

    pmu_frame estimated_frame;
    double mid_window_fracsec = 0;

    // Estimating frame
    pmu_estimate(input_signal_window, mid_window_fracsec ,estimated_frame);

    // Deinitializing pmu estimator
    pmu_deinit();

    free(input_signal_window);

    return 0;
}
