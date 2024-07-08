/*==============================================================================
  @file pmu_estimator.h

  Pmu estimator header file.

  Authors: Chemseddine Allioua, Brahim Mazighi

  Copyright (c) 2023.
  All Rights Reserved.
  Confidential and Proprietary - University of Bologna.

==============================================================================*/

#ifndef PMU_ESTIMATOR_H
#define PMU_ESTIMATOR_H

// #define NUM_CHANLS 1

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include "func_stubs.h"

#define CONFIG_FROM_INI 1
#define CONFIG_FROM_STRUCT 0

    typedef struct
    {
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

    typedef struct
    {
        float_p amp;
        float_p ph;
        float_p freq;
    } phasor;

    typedef struct
    {
        phasor synchrophasor;
        float_p rocof;
    } pmu_frame;

    // Synchrophasor Estimation Parameters
    typedef struct
    {
        uint_p win_len;
        uint_p n_cycles;
        uint_p f0;
        uint_p frame_rate;
        uint_p fs;
        uint_p n_bins;
        uint_p P;
        uint_p Q;
        bool_p iter_eipdft_enabled;
        float_p interf_trig;
        float_p df;
        float_p norm_factor;
        phasor phasor;
    } SynchrophasorEstimatorParams;

    // Dynamically allocated arrays
    typedef struct
    {
        float_p complex_p *Xf;
        float_p complex_p *Xi;
        float_p complex_p *dftbins;
        float_p *hann_window;
        float_p *signal_windows[NUM_CHANLS];
    } InternalBuffers;

    // Hanning Fourier Transform calculation constants
    typedef struct
    {
        float_p C0;
        float_p C1;
        float_p C2;
        float_p C3;
        float_p C4;
        float_p complex_p C5;
        float_p complex_p C6;
        float_p inv_norm_factor;
    } HanningTransformConstants;

    // ROCOF estimation variables
    typedef struct
    {
        float_p freq_old[NUM_CHANLS];
        float_p thresholds[3];
        float_p low_pass_coeff[3];
        float_p delay_line[NUM_CHANLS][2];
        bool_p state[NUM_CHANLS];
    } RocoFEstimationStates;

    typedef struct pmu_context
    {
        SynchrophasorEstimatorParams synch_params;
        InternalBuffers buff_params;
        HanningTransformConstants hann_params;
        RocoFEstimationStates rocof_params;
        bool_p pmu_initialized;
    } pmu_context;

    // pmu estimator initialization
    int pmu_init(pmu_context *ctx, void *cfg, bool_p config_from_ini);

    // synchrophasor, frequency, rocof estimation
    int pmu_estimate(pmu_context *ctx, float_p *in_signal_windows, float_p mid_fracsec, pmu_frame *out_frame);

    // pmu estimator deinitialization
    int pmu_deinit(pmu_context *ctx);

    // dumps quantities of a pmu_frame struct to the specified output stream
    int pmu_dump_frame(pmu_frame *frame, FILE *stream);

#ifdef __cplusplus
}
#endif

#endif /* PMU_ESTIMATOR_H */