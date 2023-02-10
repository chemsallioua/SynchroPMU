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
    unsigned int n_chanls;
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

//angle [-pi; pi] wrap in radiants
double wrap_angle(double rad_angle);

//pmu estimator initialization
int pmu_init(void* cfg);

//synchrophasor estimation
int pmu_estimate(double* in_signal_windows[], synchrophasor* out_phasor);

//pmu estimator deinitialization
int pmu_deinit();

#endif /* ITER_E_IPDFT_IMP_H */