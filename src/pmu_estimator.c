/*==============================================================================
  @file pmu_estimator.c

  Source for the implementation of the pmu estimator based on the Iterative
  Enhanced interpolated DFT Algorithm.

  Authors: Chemseddine Allioua, Brahim Mazighi

  Copyright (c) 2023.
  All Rights Reserved.
  Confidential and Proprietary - University of Bologna.

==============================================================================*/

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "iniparser.h"
#include "pmu_estimator.h"

/*CONSTANTS ==================*/
#ifndef NUM_CHANLS
#define NUM_CHANLS 1
#endif

/*LOGGING LEVEL ==============*/
#define DEBUG 3
#define INFO 2
#define ERROR 1

// change logging level here (DEBUG, INFO, ERROR)
#ifndef LOGGING_LEVEL
#define LOGGING_LEVEL ERROR
#endif

#if LOGGING_LEVEL >= ERROR
#define error(...) fprintf(stderr, __VA_ARGS__)
#else
#define error(...)
#endif
#if LOGGING_LEVEL >= INFO
#define info(...) fprintf(stdout, __VA_ARGS__)
#else
#define info(...)
#endif
#if LOGGING_LEVEL >= DEBUG
#define debug(...) fprintf(stdout, __VA_ARGS__)
#define debug_bins(...) print_bins(__VA_ARGS__)
#else
#define debug(...)
#define debug_bins(...)
#endif

/*MACROS ==================*/

#define wrap_angle(rad) (((float_p)(rad) - 2 * M_PI_p * rint((float_p)(rad) / (2 * M_PI_p)))) // wraps angle (rad) in range [-pi; pi]
#define sineasin(x) (pmue_sin(M_PI_p * (x)) / pmue_sin(ctx->hann_params.C4 * M_PI_p * (x)))
#define whDFT(x) (pmue_cexp(-I * M_PI_p * ctx->hann_params.C3 * (x)) * (ctx->hann_params.C5 * sineasin((x) - 1) + ctx->hann_params.C1 * sineasin((x)) + ctx->hann_params.C6 * sineasin((x) + 1)))
#define wf(k, f, in_phsr) ((in_phsr) * whDFT((k) - ((f) / ctx->synch_params.df)) * (ctx->hann_params.inv_norm_factor))

/*STATIC PROTOTYPES ====================*/

// The DFT implementations of a real sampled signal
static int dft_r(float_p *in_ptr, float_p complex_p *out_ptr, uint_p out_len, uint_p n_bins);
static int fft(float_p *in_ptr, float_p complex_p *out_ptr, uint_p out_len);

// helps inizializing the hann coefficients
static float_p hann(float_p *out_ptr, uint_p out_len);

// phasor and frequency estimation main functions
static void pureTone(pmu_context *ctx, float_p complex_p *Xpure, phasor phasor);
static int ipDFT(pmu_context *ctx, float_p complex_p *Xdft, phasor *phasor);
static void e_ipDFT(pmu_context *ctx, float_p complex_p *Xdft, phasor *out_phasor);
static void iter_e_ipDFT(pmu_context *ctx, float_p complex_p *dftbins, float_p complex_p *Xi, float_p complex_p *Xf, phasor *f_phsr);

// phasor and frequency estimation helping functions
inline static void find3LargestIndx(float_p arr[], int size, uint_p *km, uint_p *kl, uint_p *kr);

// pmu estimator configuration functions
static int config_estimator(pmu_context *ctx, void *cfg, bool_p config_from_ini);
static int check_config_validity(pmu_context *ctx);
static int config_from_file(pmu_context *ctx, char *ini_file_name);

// prints the bins and their index and frequency
static void print_bins(float_p complex_p *bins, int n_bins, float_p df, char *str);

/*IMPLEMENTATION ====================*/

// pmu initialization function implementation
int pmu_init(pmu_context *ctx, void *cfg, bool_p config_from_ini)
{

    // check if the pmu estimator is already initialized
    if (ctx->pmu_initialized)
    {
        error("[%s] ERROR: pmu estimator already initialized\n", __FUNCTION__);
        return -1;
    }
    if (config_estimator(ctx, cfg, config_from_ini))
    {
        error("[%s] ERROR: pmu estimator configuration failed\n", __FUNCTION__);
        return -1;
    }

    info("[%s] Initializing pmu estimator\n", __FUNCTION__);

    // allocate memory for the global arrays
    if (NULL == (ctx->buff_params.Xf = malloc(ctx->synch_params.n_bins * sizeof(float_p complex_p))))
    {
        error("[%s] ERROR: ctx->buff_params.Xf memory allocation failed\n", __FUNCTION__);
        return -1;
    }
    if (NULL == (ctx->buff_params.Xi = malloc(ctx->synch_params.n_bins * sizeof(float_p complex_p))))
    {
        error("[%s] ERROR: ctx->buff_params.Xi memory allocation failed\n", __FUNCTION__);
        return -1;
    }
    if (NULL == (ctx->buff_params.dftbins = malloc(ctx->synch_params.win_len * sizeof(float_p complex_p))))
    {
        error("[%s] ERROR: ctx->buff_params.dftbins memory allocation failed\n", __FUNCTION__);
        return -1;
    }
    if (NULL == (ctx->buff_params.hann_window = malloc(ctx->synch_params.win_len * sizeof(float_p))))
    {
        error("[%s] ERROR: ctx->buff_params.hann_window memory allocation failed\n", __FUNCTION__);
        return -1;
    }
    uint_p i;
    for (i = 0; i < NUM_CHANLS; i++)
    {
        if (NULL == (ctx->buff_params.signal_windows[i] = (float_p *)malloc(ctx->synch_params.win_len * sizeof(float_p))))
        {
            error("[%s] ERROR: ctx->buff_params.signal_windows memory allocation failed\n", __FUNCTION__);
            return -1;
        }
    }

    // initialize the global arrays
    for (i = 0; i < NUM_CHANLS; i++)
    {
        ctx->rocof_params.delay_line[i][0] = 0.0;
        ctx->rocof_params.delay_line[i][1] = 0.0;
        ctx->rocof_params.freq_old[i] = ctx->synch_params.f0;
        ctx->rocof_params.state[i] = 0;
    }

    // initialize the hann coefficients
    ctx->synch_params.norm_factor = hann(ctx->buff_params.hann_window, ctx->synch_params.win_len);

    // initialize Hanning FT calculation constants
    ctx->hann_params.C0 = -0.25;
    ctx->hann_params.C1 = 0.5;
    ctx->hann_params.C2 = -0.25;
    ctx->hann_params.C3 = ((float_p)ctx->synch_params.win_len - 1) / (float_p)ctx->synch_params.win_len;
    ctx->hann_params.C4 = 1.0 / (float_p)ctx->synch_params.win_len;
    ctx->hann_params.C5 = ctx->hann_params.C0 * pmue_cexp(I * M_PI * ctx->hann_params.C3);
    ctx->hann_params.C6 = ctx->hann_params.C2 * pmue_cexp(-I * M_PI * ctx->hann_params.C3);
    ctx->hann_params.inv_norm_factor = 1.0 / ctx->synch_params.norm_factor;

    // pmu estimator is initialized successfully
    ctx->pmu_initialized = 1;

    info("[%s] Pmu estimator initialized successfully\n", __FUNCTION__);
    return 0;
}

// pmu estimation function implementation
int pmu_estimate(pmu_context *ctx, float_p *in_signal_windows, float_p mid_fracsec, pmu_frame *out_frame)
{

    info("[%s] pmu_estimate() started\n", __FUNCTION__);

    // check if pmu estimator is initialized
    if (!ctx->pmu_initialized)
    {
        error("[%s] ERROR: pmu estimator not initialized, first initialize with pmu_init function\n", __FUNCTION__);
        return -1;
    }

    // input signal windowing
    uint_p j, chnl;
    for (chnl = 0; chnl < NUM_CHANLS; chnl++)
    {
        for (j = 0; j < ctx->synch_params.win_len; j++)
        {
            ctx->buff_params.signal_windows[chnl][j] = (*((in_signal_windows + chnl * ctx->synch_params.win_len) + j)) * ctx->buff_params.hann_window[j];
            // info("in_signal_windows[%u][%u]§; %f\n", chnl, j, *((in_signal_windows + chnl * ctx->synch_params.win_len) + j));
        }
    }

    info("[%s] windowing on all channels done successfully\n", __FUNCTION__);

    // synchrophasor estimation for each channel
    for (chnl = 0; chnl < NUM_CHANLS; chnl++)
    {

        // compuute DFT of input signal
        pmue_fft_r((float_p *)ctx->buff_params.signal_windows[chnl], ctx->buff_params.dftbins, ctx->synch_params.win_len, ctx->synch_params.n_bins);

        debug_bins(ctx->buff_params.dftbins, ctx->synch_params.n_bins, ctx->synch_params.df, "Input Signal DFT BINS");

        // perform enhanced interpolated DFT
        e_ipDFT(ctx, ctx->buff_params.dftbins, &ctx->synch_params.phasor);
        pureTone(ctx, ctx->buff_params.Xf, ctx->synch_params.phasor);

        float_p E_diff = 0;
        float_p E = 0;

        // compute energy of both input signal and interference frequencies
        for (j = 0; j < ctx->synch_params.n_bins; j++)
        {

            ctx->buff_params.Xi[j] = ctx->buff_params.dftbins[j] - ctx->buff_params.Xf[j];

            if (ctx->synch_params.iter_eipdft_enabled)
            {
                E_diff += pmue_cabs(ctx->buff_params.Xi[j] * ctx->buff_params.Xi[j]);
                E += pmue_cabs(ctx->buff_params.dftbins[j] * ctx->buff_params.dftbins[j]);
            }
        }

        debug("Energy of Signal Spectrum (multiplied by epsilon) = %lf | Energy of Difference= %lf\n", E * ctx->synch_params.interf_trig, E_diff);

        // check if interference is present to trigger iterative enhanced interpolated DFT
        if (ctx->synch_params.iter_eipdft_enabled)
        {
            if (E_diff > ctx->synch_params.interf_trig * E)
            {
                iter_e_ipDFT(ctx, ctx->buff_params.dftbins, ctx->buff_params.Xi, ctx->buff_params.Xf, &ctx->synch_params.phasor);
            }
        }

        // two state rocof estimation
        float_p rocof = (ctx->synch_params.phasor.freq - ctx->rocof_params.freq_old[chnl]) * (float_p)ctx->synch_params.frame_rate;
        float_p rocof_der = (rocof - ctx->rocof_params.delay_line[chnl][0]) * (float_p)ctx->synch_params.frame_rate;

        // update state
        ctx->rocof_params.state[chnl] = (!ctx->rocof_params.state[chnl] && (pmue_fabs(rocof) > ctx->rocof_params.thresholds[0] || pmue_fabs(rocof_der) > ctx->rocof_params.thresholds[1])) ? 1 : ctx->rocof_params.state[chnl];
        ctx->rocof_params.state[chnl] = (ctx->rocof_params.state[chnl] && pmue_fabs(rocof) < ctx->rocof_params.thresholds[2]) ? 0 : ctx->rocof_params.state[chnl];

        // apply low pass filter if state is 0
        if (!ctx->rocof_params.state[chnl])
        {
            out_frame[chnl].rocof = ctx->rocof_params.low_pass_coeff[1] * rocof +
                                    ctx->rocof_params.low_pass_coeff[2] * ctx->rocof_params.delay_line[chnl][0] -
                                    ctx->rocof_params.low_pass_coeff[0] * ctx->rocof_params.delay_line[chnl][1];
        }
        else
        {
            out_frame[chnl].rocof = rocof;
        }

        // update old frequency
        ctx->rocof_params.freq_old[chnl] = ctx->synch_params.phasor.freq;

        // update delay line
        ctx->rocof_params.delay_line[chnl][0] = rocof;                 // x(n-1)
        ctx->rocof_params.delay_line[chnl][1] = out_frame[chnl].rocof; // y(n-1)

        // populate output frame
        out_frame[chnl].synchrophasor.freq = ctx->synch_params.phasor.freq;
        out_frame[chnl].synchrophasor.amp = 2 * ctx->synch_params.phasor.amp / ctx->synch_params.norm_factor;
        out_frame[chnl].synchrophasor.ph = wrap_angle(ctx->synch_params.phasor.ph - 2 * M_PI_p * ctx->synch_params.f0 * mid_fracsec);

        info("[%s] Estimated on channel: %d\n", __FUNCTION__, chnl);
    }

    return 0;
}

// pmu deinitialization fuunction implementation
int pmu_deinit(pmu_context *ctx)
{

    info("[%s] Deinitializing pmu estimator\n", __FUNCTION__);

    // check if pmu estimator is initialized
    if (!ctx->pmu_initialized)
    {
        error("[%s] ERROR: pmu estimator not initialized, first initialize with pmu_init function\n", __FUNCTION__);
        return -1;
    }

    // free all allocated memory
    free(ctx->buff_params.Xf);
    free(ctx->buff_params.Xi);
    free(ctx->buff_params.dftbins);
    free(ctx->buff_params.hann_window);

    uint_p i;
    for (i = 0; i < NUM_CHANLS; i++)
    {
        free(ctx->buff_params.signal_windows[i]);
    }

    info("[%s] Pmu estimator deinitialized successfully\n", __FUNCTION__);

    // set pmu estimator to uninitialized state
    ctx->pmu_initialized = 0;
    return 0;
}

// dumps quantities of a pmu_frame struct to the specified output stream
int pmu_dump_frame(pmu_frame *frame, FILE *stream)
{

    if (frame == NULL || stream == NULL)
    {
        error("Error: NULL pointer passed as argument\n");
        return -1;
    }

    int written = fprintf(stream, "[Synchrophasor] amplitude: %lf, phase: %lf, frequency: %lf, rocof: %lf\n",
                          frame->synchrophasor.amp, wrap_angle(frame->synchrophasor.ph), frame->synchrophasor.freq, frame->rocof);
    if (written < 0)
    {
        error("Error: failed to write to stream\n");
        return -1;
    }

    return 0;
}

// pmu estimator configuration functions
static int config_estimator(pmu_context *ctx, void *cfg, bool_p config_from_ini)
{

    info("[%s] Configurating pmu estimator\n", __FUNCTION__);

    // pmu estimator parameters initialization from input config
    if (config_from_ini)
    {
        char *ini_file_name = (char *)cfg;
        if (config_from_file(ctx, ini_file_name))
        {
            error("[%s] ERROR: pmu estimator configuration failed\n", __FUNCTION__);
            return -1;
        }
    }
    else
    {
        estimator_config *config = (estimator_config *)cfg;

        ctx->synch_params.fs = config->fs;
        ctx->synch_params.f0 = config->f0;
        ctx->synch_params.n_cycles = config->n_cycles;
        ctx->synch_params.frame_rate = config->frame_rate;
        ctx->synch_params.n_bins = config->n_bins;
        ctx->synch_params.P = config->P;
        ctx->synch_params.Q = config->Q;
        ctx->synch_params.interf_trig = config->interf_trig;
        ctx->synch_params.iter_eipdft_enabled = config->iter_eipdft;

        ctx->rocof_params.thresholds[0] = config->rocof_thresh[0];
        ctx->rocof_params.thresholds[1] = config->rocof_thresh[1];
        ctx->rocof_params.thresholds[2] = config->rocof_thresh[2];

        ctx->rocof_params.low_pass_coeff[0] = config->rocof_low_pass_coeffs[0];
        ctx->rocof_params.low_pass_coeff[1] = config->rocof_low_pass_coeffs[1];
        ctx->rocof_params.low_pass_coeff[2] = config->rocof_low_pass_coeffs[2];
    }

    ctx->synch_params.win_len = ctx->synch_params.n_cycles * ctx->synch_params.fs / ctx->synch_params.f0;
    ctx->synch_params.df = (float_p)ctx->synch_params.fs / (float_p)ctx->synch_params.win_len;

    if (check_config_validity(ctx))
    {
        error("[%s] ERROR: pmu estimator configuration failed, config values not valid\n", __FUNCTION__);
        return -1;
    }

    debug("\n[%s] Configuration: win_len: %u, n_cycles: %u, fs: %u, f0: %u,\
            frame_rate: %u\n n_bins: %u, P: %u, Q: %u, interf_trig: %f\n df: %f,\
            thresholds: %f, %f, %f\n low_pass_coeff: %f, %f, %f\n\n",
          __FUNCTION__, win_len, ctx->synch_params.n_cycles, ctx->synch_params.fs, ctx->synch_params.f0, ctx->synch_params.frame_rate, ctx->synch_params.n_bins, ctx->synch_params.P, ctx->synch_params.Q,
          ctx->synch_params.interf_trig, ctx->synch_params.df, ctx->rocof_params.thresholds[0], ctx->rocof_params.thresholds[1], ctx->rocof_params.thresholds[2],
          ctx->rocof_params.low_pass_coeff[0], ctx->rocof_params.low_pass_coeff[1], ctx->rocof_params.low_pass_coeff[2]);

    info("[%s] Configuration done successfully \n", __FUNCTION__);

    return 0;
}

static int check_config_validity(pmu_context *ctx)
{

    if (ctx->synch_params.f0 != 50 && ctx->synch_params.f0 != 60)
    {
        error("[%s] ERROR: nominal frequency: %u not correctly set, (allowed values 50 or 60)\n", __FUNCTION__, ctx->synch_params.f0);
        return -1;
    }

    if (ctx->synch_params.n_cycles <= 0 || ctx->synch_params.n_cycles > 50)
    {
        error("[%s] ERROR: window length: %u is not correctly set, must be non-zero, positive and smaller than 50\n", __FUNCTION__, ctx->synch_params.n_cycles);
        return -1;
    }
    if (ctx->synch_params.fs <= 0 || fmod(ctx->synch_params.fs, ctx->synch_params.f0) != 0.0)
    {
        error("[%s] ERROR: sample rate: %u is not correctly set, fs must be non-zero, positive and can be devided by f0: %u\n", __FUNCTION__, ctx->synch_params.fs, ctx->synch_params.f0);
        return -1;
    }

    if (ctx->synch_params.frame_rate <= 0)
    {
        error("[%s] ERROR: frame rate: %u is not correctly set, must be non-zero and positive\n", __FUNCTION__, ctx->synch_params.frame_rate);
        return -1;
    }
    if (ctx->synch_params.n_bins < ctx->synch_params.n_cycles + 2 || ctx->synch_params.n_bins > (ctx->synch_params.n_cycles * ctx->synch_params.fs / ctx->synch_params.f0) / 2)
    {
        error("[%s] ERROR: number of dft bins: %u for estimation is not correctly set, must be greater than or equal (number of cycles +2): %d and smaller than or equal half the window length: %d \n", __FUNCTION__, ctx->synch_params.n_bins, ctx->synch_params.n_cycles + 2, (ctx->synch_params.n_cycles * ctx->synch_params.fs / ctx->synch_params.f0) / 2);
        return -1;
    }
    if (ctx->synch_params.P <= 0)
    {
        error("[%s] ERROR: ipdft iterations: %u is not correctly set, must be non-zero and positive\n", __FUNCTION__, ctx->synch_params.P);
        return -1;
    }
    if (ctx->synch_params.interf_trig <= 0 || ctx->synch_params.interf_trig > 1)
    {
        error("[%s] ERROR: interference threshold: %lf is not correctly set, must be between ]0,1]\n", __FUNCTION__, ctx->synch_params.interf_trig);
        return -1;
    }
    if (NUM_CHANLS < 1 || fmod(NUM_CHANLS, 1) != 0)
    {
        error("[%s] ERROR: number of channels: %u is not correctly set, must be an integer and at least 1\n", __FUNCTION__, NUM_CHANLS);
        return -1;
    }

    if (ctx->rocof_params.thresholds[0] <= 0)
    {
        error("[%s] ERROR: rocof threshold 1: %lf is not set, must be non-zero and positive\n", __FUNCTION__, ctx->rocof_params.thresholds[0]);
        return -1;
    }
    if (ctx->rocof_params.thresholds[1] <= 0)
    {
        error("[%s] ERROR: rocof threshold 2: %lf is not set, must be non-zero and positive\n", __FUNCTION__, ctx->rocof_params.thresholds[1]);
        return -1;
    }
    if (ctx->rocof_params.thresholds[2] <= 0)
    {
        error("[%s] ERROR: rocof threshold 3: %lf is not set, must be non-zero and positive\n", __FUNCTION__, ctx->rocof_params.thresholds[2]);
        return -1;
    }

    return 0;
}

static int config_from_file(pmu_context *ctx, char *ini_file_name)
{

    dictionary *ini;

    ini = iniparser_load(ini_file_name);
    if (ini == NULL)
    {
        error("[%s] Error: cannot parse file: %s\n", __FUNCTION__, ini_file_name);
        return -1;
    }

    ctx->synch_params.n_cycles = (int)iniparser_getint(ini, "signal:n_cycles", 0);
    ctx->synch_params.fs = (int)iniparser_getint(ini, "signal:sample_rate", 0);
    ctx->synch_params.f0 = (int)iniparser_getint(ini, "signal:nominal_freq", 0);
    ctx->synch_params.frame_rate = (int)iniparser_getint(ini, "synchrophasor:frame_rate", 0);
    ctx->synch_params.n_bins = (int)iniparser_getint(ini, "synchrophasor:number_of_dft_bins", 0);
    ctx->synch_params.P = (int)iniparser_getint(ini, "synchrophasor:ipdft_iterations", 0);
    ctx->synch_params.Q = (int)iniparser_getint(ini, "synchrophasor:iter_e_ipdft_iterations", 0);
    ctx->synch_params.interf_trig = (float_p)iniparser_getdouble(ini, "synchrophasor:interference_threshold", 0);
    ctx->synch_params.iter_eipdft_enabled = (bool_p)iniparser_getboolean(ini, "synchrophasor:iter_e_ipdft_enable", 0);

    ctx->rocof_params.thresholds[0] = (float_p)iniparser_getdouble(ini, "rocof:threshold_1", 0);
    ctx->rocof_params.thresholds[1] = (float_p)iniparser_getdouble(ini, "rocof:threshold_2", 0);
    ctx->rocof_params.thresholds[2] = (float_p)iniparser_getdouble(ini, "rocof:threshold_3", 0);

    ctx->rocof_params.low_pass_coeff[0] = (float_p)iniparser_getdouble(ini, "rocof:low_pass_filter_1", 0);
    ctx->rocof_params.low_pass_coeff[1] = (float_p)iniparser_getdouble(ini, "rocof:low_pass_filter_2", 0);
    ctx->rocof_params.low_pass_coeff[2] = (float_p)iniparser_getdouble(ini, "rocof:low_pass_filter_3", 0);

    iniparser_freedict(ini);

    return 0;
}

// performs the DFT of a real sampled signal
static int dft_r(float_p *in_ptr, float_p complex_p *out_ptr, uint_p out_len, uint_p n_bins)
{

    if (((out_len & (out_len - 1)) == 0))
    {
        float_p g_fft_in[out_len];
        memcpy(g_fft_in, in_ptr, out_len * sizeof(float_p));
        fft(g_fft_in, out_ptr, out_len);
    }
    else
    {

        uint_p k, n;

        for (k = 0; k < n_bins; ++k)
        {
            out_ptr[k] = 0;
            for (n = 0; n < out_len; ++n)
                out_ptr[k] += (in_ptr[n] * pmue_cexp(-I * ((n * k * M_PI_p * 2 / (float_p)out_len))));
        }
    }

    return 0;
}
static int fft(float_p *in_ptr, float_p complex_p *out_ptr, uint_p out_len)
{

    uint_p half_len = out_len / 2;
    uint_p i, k;

    if (out_len == 1)
    {
        *out_ptr = in_ptr[0];
    }
    else
    {
        float_p even[half_len], odd[half_len];
        float_p complex_p even_out[half_len], odd_out[half_len];

        // Split input into even and odd indexed elements
        for (i = 0; i < half_len; i++)
        {
            even[i] = in_ptr[2 * i];
            odd[i] = in_ptr[2 * i + 1];
        }

        // Compute FFT of even and odd indexed elements recursively
        fft(even, even_out, half_len);
        fft(odd, odd_out, half_len);

        // Combine results from even and odd indexed elements
        for (k = 0; k < half_len; k++)
        {
            float_p complex_p t = pmue_cexp(-I * 2 * M_PI * k / out_len) * odd_out[k];
            out_ptr[k] = even_out[k] + t;
            out_ptr[k + half_len] = even_out[k] - t;
        }
    }
    return 0;
}

// helps inizializing the hann coefficients
static float_p hann(float_p *out_ptr, uint_p out_len)
{

    float_p norm_fact = 0;
    uint_p i = 0;
    for (i = 0; i < out_len; i++)
    {
        out_ptr[i] = 0.5 * (1 - pmue_cos(2 * M_PI_p * i / out_len));
        norm_fact += out_ptr[i];
    }
    return norm_fact;
}

// phasor and frequency estimation main functions
static void pureTone(pmu_context *ctx, float_p complex_p *Xpure, phasor phasor)
{
    debug("\n[pureTone] ===============================================\n");
    uint_p i;
    float_p complex_p phsr_0 = phasor.amp * pmue_cexp(I * phasor.ph);
    float_p complex_p phsr_1 = phasor.amp * pmue_cexp(-I * phasor.ph);
    for (i = 0; i < ctx->synch_params.n_bins; i++)
    {
        Xpure[i] = wf(i, phasor.freq, phsr_0) + wf(i, -phasor.freq, phsr_1);
    }

    debug("freq: %0.3lf | ampl: %0.3lf | phse: %0.3lf\n", phasor.freq, phasor.amp, phasor.ph);
    debug_bins(Xpure, ctx->synch_params.n_bins, ctx->synch_params.df, "DFT BINS PURE TONE");
    debug("[END pureTone] ================================================\n\n");
}

static int ipDFT(pmu_context *ctx, float_p complex_p *Xdft, phasor *phasor)
{

    uint_p j, k1, k2, k3;
    float_p Xdft_mag[ctx->synch_params.n_bins]; // magnitude of dft

    debug("\n[ipDFT] ===============================================\n");

    debug_bins(Xdft, ctx->synch_params.n_bins, ctx->synch_params.df, "DFT BINS IN ipDFT");

    for (j = 0; j < ctx->synch_params.n_bins; j++)
    {
        Xdft_mag[j] = pmue_cabs(Xdft[j]);
    }

    find3LargestIndx(Xdft_mag, ctx->synch_params.n_bins, &k1, &k2, &k3);

    debug("[%s] k1: %d, k2: %d, k3: %d\n", __FUNCTION__, k1, k2, k3);

    float_p delta_corr = 2 * (Xdft_mag[k3] - Xdft_mag[k2]) / (Xdft_mag[k2] + Xdft_mag[k3] + 2 * Xdft_mag[k1]);

    debug("[%s] delta_corr: %lf\n", __FUNCTION__, delta_corr);
    phasor->freq = (k1 + delta_corr) * ctx->synch_params.df;

    if (pmue_fabs(delta_corr) <= 10e-12)
    {

        phasor->amp = Xdft_mag[k1];
        phasor->ph = pmue_carg(Xdft[k1]);

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n", __FUNCTION__, phasor->freq, phasor->amp, phasor->ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 1;
    }
    else
    {
        phasor->amp = Xdft_mag[k1] * pmue_fabs((delta_corr * delta_corr - 1) * (M_PI_p * delta_corr) / pmue_sin(M_PI_p * delta_corr));
        phasor->ph = pmue_carg(Xdft[k1]) - M_PI_p * delta_corr;

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n", __FUNCTION__, phasor->freq, phasor->amp, phasor->ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 0;
    }
}

static void e_ipDFT(pmu_context *ctx, float_p complex_p *Xdft, phasor *out_phasor)
{

    debug("\n[e_ipDFT] ===============================================\n");

    phasor phsr = *out_phasor;

    if (!ipDFT(ctx, Xdft, &phsr))
    {
        uint_p j, p;

        float_p complex_p X_neg;
        float_p complex_p X_pos[ctx->synch_params.n_bins];
        float_p complex tmp_phsr;

        for (p = 0; p < ctx->synch_params.P; p++)
        {

            debug("\n[e_ipDFT ITERATION: %d] ------------\n", p + 1);

            tmp_phsr = phsr.amp * pmue_cexp(-I * phsr.ph);

            for (j = 0; j < ctx->synch_params.n_bins; j++)
            {
                X_neg = wf(j, -phsr.freq, tmp_phsr);
                X_pos[j] = Xdft[j] - X_neg;
            }

            if (ipDFT(ctx, X_pos, &phsr))
            {
                break;
            }
            debug("\nEND e_ipDFT ITERATION --------------------------\n");
        }
    }
    debug("\n[END e_ipDFT]========================================================\n\n");

    *out_phasor = phsr;
}

static void iter_e_ipDFT(pmu_context *ctx, float_p complex_p *dftbins, float_p complex_p *Xi, float_p complex_p *Xf, phasor *f_phsr)
{

    phasor f_phasor = *f_phsr;
    phasor i_phasor;
    float_p complex_p Xi_pure[ctx->synch_params.n_bins];

    debug("\n[iter-e-ipDFT] ###############################################\n");
    uint_p i, j;
    for (i = 0; i < ctx->synch_params.Q; i++)
    {
        debug("\n[iter-e-ipDFT ITERATION: %d] ------------\n", i + 1);

        e_ipDFT(ctx, Xi, &i_phasor);
        pureTone(ctx, Xi_pure, i_phasor);
        for (j = 0; j < ctx->synch_params.n_bins; j++)
        {
            Xf[j] = dftbins[j] - Xi_pure[j];
        }
        e_ipDFT(ctx, Xf, &f_phasor);

        if (i < ctx->synch_params.Q - 1)
        {

            pureTone(ctx, Xf, f_phasor);

            for (j = 0; j < ctx->synch_params.n_bins; j++)
            {
                Xi[j] = dftbins[j] - Xf[j];
            }
        }

        debug("\nEND iter-e-ipDFT ITERATION --------------------------\n");
    }
    debug("\n[END iter-e-ipDFT] ##############################################\n\n");

    *f_phsr = f_phasor;
}

inline static void find3LargestIndx(float_p arr[], int size, uint_p *km, uint_p *kl, uint_p *kr)
{
    int max_indx = 0;

    for (int i = 1; i < size; i++)
    {
        if (arr[i] > arr[max_indx])
        {
            max_indx = i;
        }
    }

    *km = max_indx;
    *kl = (max_indx == 0) ? 0 : max_indx - 1;
    *kr = (max_indx == size - 1) ? max_indx - 1 : max_indx + 1;
}

// prints the bins and their index and frequency
static void print_bins(float_p complex_p *bins, int n_bins, float_p df, char *str)
{

    debug("\n--%s---------------  ---  --  -\n Indx", str);
    for (int i = 0; i < n_bins; i++)
    {
        debug("|%6d", i);
    }
    debug("|\n Freq");
    for (int i = 0; i < n_bins; i++)
    {
        debug("|%6.1f", i * df);
    }
    debug("|\n Bins");
    for (int i = 0; i < n_bins; i++)
    {
        debug("|%6.1f", pmue_cabs(bins[i]));
    }
    debug("|\n---------------------------  ---  --  -\n\n");
}
