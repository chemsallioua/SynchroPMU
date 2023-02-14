/*==============================================================================
  @file pmu_estimator.h

  Pmu estimator header file. 

  Authors: Chemseddine Allioua, Brahim Mazighi 

  Copyright (c) 2023.
  All Rights Reserved.
  Confidential and Proprietary - University of Bologna.

==============================================================================*/ 
#include <stdio.h>

#ifndef PMU_ESTIMATOR_H
#define PMU_ESTIMATOR_H

#ifdef __cplusplus
extern "C" {
#endif

//#define DEBUG 0

#if DEBUG
#define debug(...) printf(__VA_ARGS__)
#define debug_bins(...) print_bins(__VA_ARGS__)
#else
#define debug(...)
#define debug_bins(...)
#endif

#define CONFIG_FROM_INI 1
#define CONFIG_FROM_STRUCT 0

typedef struct {
    unsigned int n_cycles;
    unsigned int f0;
    unsigned int frame_rate;
    unsigned int fs;
    unsigned int n_bins;
    unsigned int P;
    unsigned int Q;
    _Bool iter_eipdft;
    double interf_trig;
    double rocof_thresh[3];
    double rocof_low_pass_coeffs[3];
} estimator_config;

typedef struct{
    double amp;
    double ph;
    double freq;
}phasor;

typedef struct{
    phasor synchrophasor;
    double rocof;
}pmu_frame;

//pmu estimator initialization
int pmu_init(void* cfg, _Bool config_from_ini);

//synchrophasor, frequency, rocof estimation
int pmu_estimate(double *in_signal_windows, double mid_fracsec ,pmu_frame* out_frame);

//pmu estimator deinitialization
int pmu_deinit();

// dumps quantities of a pmu_frame struct to the specified output stream
int pmu_dump_frame(pmu_frame *frame, FILE *stream);

#ifdef __cplusplus
}
#endif

#endif /* ITER_E_IPDFT_IMP_H */