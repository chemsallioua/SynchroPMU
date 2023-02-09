/*==============================================================================
  @file iter_e_ipdft_imp.h

  Pmu estimator header file. 

  Copyright (c) 2023.
  All Rights Reserved.
  Confidential and Proprietary - University of Bologna.

==============================================================================*/ 

#ifndef ITER_E_IPDFT_IMP_H
#define ITER_E_IPDFT_IMP_H

#include <stdio.h>
#include <complex.h>
#include <math.h>

#define DEBUG 0

#if DEBUG
#define debug(...) printf(__VA_ARGS__)
#define debug_bins(...) print_bins(__VA_ARGS__)
#else
#define debug(...)
#define debug_bins(...)
#endif

typedef struct {
    unsigned int win_len;
    double fs;
    unsigned int n_bins;
    unsigned int P;
    unsigned int Q;
    double interf_trig;
} estimator_config;

typedef struct{
    double amp;
    double ph;
    double freq;
}synchrophasor;

inline float wrap_angle(float rad_angle){
	float temp = fmod(rad_angle + M_PI, 2*M_PI);
	if(temp < 0.0){
        temp += 2.0*M_PI;
		}

	return temp - M_PI;
}

//pmu estimator initialization
int pmu_init(void* cfg);

//synchrophasor estimation
int pmu_estimate(double * signal_window, synchrophasor* out_phasor);

//pmu estimator deinitialization
int pmu_deinit();

#endif /* ITER_E_IPDFT_IMP_H */