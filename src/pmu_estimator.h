/*==============================================================================
  @file pmu_estimator.h

  Pmu estimator header file. 

  Authors: Chemseddine Allioua, Brahim Mazighi 

  Copyright (c) 2023.
  All Rights Reserved.
  Confidential and Proprietary - University of Bologna.

==============================================================================*/ 
#include <stdio.h>
#include "func_stubs.h"

#ifndef PMU_ESTIMATOR_H
#define PMU_ESTIMATOR_H

#ifdef __cplusplus
extern "C" {
#endif

//#define DEBUG 1

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
    uint_p n_cycles;
    uint_p f0;
    uint_p frame_rate;
    uint_p fs;
    uint_p n_bins;
    uint_p P;
    uint_p Q;
    bool_p iter_eipdft;
    float_p interf_trig;
    float_p rocof_thresh[3];
    float_p rocof_low_pass_coeffs[3];
} estimator_config;

typedef struct{
    float_p amp;
    float_p ph;
    float_p freq;
}phasor;

typedef struct{
    phasor synchrophasor;
    float_p rocof;
}pmu_frame;

//pmu estimator initialization
int pmu_init(void* cfg, bool_p config_from_ini);

//synchrophasor, frequency, rocof estimation
int pmu_estimate(float_p *in_signal_windows, float_p mid_fracsec ,pmu_frame* out_frame);

//pmu estimator deinitialization
int pmu_deinit();

// dumps quantities of a pmu_frame struct to the specified output stream
int pmu_dump_frame(pmu_frame *frame, FILE *stream);

#ifdef __cplusplus
}
#endif

#endif /* PMU_ESTIMATOR_H */