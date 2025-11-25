/**
 * @file pmu_estimator.h
 * @brief PMU (Phasor Measurement Unit) Estimator Header File
 * 
 * This header file provides the interface for the PMU estimator library, which implements
 * the Iterative Enhanced Interpolated DFT (Discrete Fourier Transform) Synchrophasor 
 * Estimation Algorithm.
 * 
 * The library provides functionality for:
 * - Synchrophasor estimation (amplitude, phase, frequency)
 * - Rate of Change of Frequency (ROCOF) estimation
 * - Multi-channel signal processing
 * - Configurable estimation parameters
 * 
 * @author Chemseddine Allioua, Brahim Mazighi
 * @copyright Copyright (c) 2023. All Rights Reserved.
 *            Confidential and Proprietary - University of Bologna.
 * @version 1.7.0
 */

#ifndef PMU_ESTIMATOR_H
#define PMU_ESTIMATOR_H

// #define NUM_CHANLS 1

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include "func_stubs.h"

/** @brief Configuration source: Load configuration from INI file */
#define CONFIG_FROM_INI 1
/** @brief Configuration source: Load configuration from structure */
#define CONFIG_FROM_STRUCT 0

    /**
     * @struct estimator_config
     * @brief Configuration structure for the PMU estimator
     * 
     * This structure contains all configuration parameters needed to initialize
     * and configure the PMU estimator. Can be loaded from an INI file or passed
     * directly as a structure.
     */
    typedef struct
    {
        uint_p n_cycles;                    /**< Number of cycles of the fundamental frequency in the observation window */
        uint_p f0;                          /**< Nominal fundamental frequency (Hz) */
        uint_p frame_rate;                  /**< Output frame rate (frames per second) */
        uint_p fs;                          /**< Sampling frequency (Hz) */
        uint_p n_bins;                      /**< Number of DFT bins to compute */
        uint_p P;                           /**< Number of iterations for the iterative algorithm */
        uint_p Q;                           /**< Number of interference tones for enhanced algorithm */
        bool_p iter_eipdft;                 /**< Enable/disable iterative enhanced interpolated DFT */
        float_p interf_trig;                /**< Interference detection trigger threshold */
        float_p rocof_thresh[3];            /**< ROCOF thresholds for different conditions */
        float_p rocof_low_pass_coeffs[3];   /**< Low-pass filter coefficients for ROCOF estimation */
    } estimator_config;

    /**
     * @struct phasor
     * @brief Represents a synchrophasor (complex voltage/current measurement)
     * 
     * A phasor is a complex number representation of a sinusoidal signal,
     * containing amplitude, phase, and frequency information.
     */
    typedef struct
    {
        float_p amp;   /**< Phasor amplitude (magnitude) */
        float_p ph;    /**< Phasor phase angle (radians) */
        float_p freq;  /**< Frequency (Hz) */
    } phasor;

    /**
     * @struct pmu_frame
     * @brief Complete PMU output frame
     * 
     * Contains the synchrophasor estimate and the Rate of Change of Frequency (ROCOF).
     * This is the primary output structure of the pmu_estimate() function.
     */
    typedef struct
    {
        phasor synchrophasor;  /**< Estimated synchrophasor */
        float_p rocof;         /**< Rate of Change of Frequency (Hz/s) */
    } pmu_frame;

    /**
     * @struct SynchrophasorEstimatorParams
     * @brief Internal parameters for synchrophasor estimation
     * 
     * This structure holds the processed configuration parameters used internally
     * by the estimation algorithm. These are derived from the estimator_config.
     */
    typedef struct
    {
        uint_p win_len;               /**< Window length in samples */
        uint_p n_cycles;              /**< Number of cycles in observation window */
        uint_p f0;                    /**< Nominal frequency (Hz) */
        uint_p frame_rate;            /**< Frame rate (frames/s) */
        uint_p fs;                    /**< Sampling frequency (Hz) */
        uint_p n_bins;                /**< Number of DFT bins */
        uint_p P;                     /**< Iteration count for iterative algorithm */
        uint_p Q;                     /**< Number of interference tones */
        bool_p iter_eipdft_enabled;   /**< Flag: iterative enhanced ipDFT enabled */
        float_p interf_trig;          /**< Interference detection threshold */
        float_p df;                   /**< Frequency resolution (Hz) */
        float_p norm_factor;          /**< Normalization factor for DFT */
        phasor phasor;                /**< Current phasor estimate */
    } SynchrophasorEstimatorParams;

    /**
     * @struct InternalBuffers
     * @brief Dynamically allocated internal buffers
     * 
     * Contains pointers to dynamically allocated arrays used during the
     * estimation process. These are allocated during initialization and
     * freed during deinitialization.
     */
    typedef struct
    {
        float_p complex_p *Xf;                      /**< DFT output for final phasor */
        float_p complex_p *Xi;                      /**< DFT output for intermediate calculations */
        float_p complex_p *dftbins;                 /**< DFT bins array */
        float_p *hann_window;                       /**< Hanning window coefficients */
        float_p *signal_windows[NUM_CHANLS];        /**< Signal window buffers for each channel */
    } InternalBuffers;

    /**
     * @struct HanningTransformConstants
     * @brief Precomputed constants for Hanning window Fourier Transform
     * 
     * These constants are calculated during initialization and used to optimize
     * the Hanning window DFT calculations during estimation.
     */
    typedef struct
    {
        float_p C0;                  /**< Constant C0 for Hanning transform */
        float_p C1;                  /**< Constant C1 for Hanning transform */
        float_p C2;                  /**< Constant C2 for Hanning transform */
        float_p C3;                  /**< Constant C3 for Hanning transform */
        float_p C4;                  /**< Constant C4 for Hanning transform */
        float_p complex_p C5;        /**< Complex constant C5 for Hanning transform */
        float_p complex_p C6;        /**< Complex constant C6 for Hanning transform */
        float_p inv_norm_factor;     /**< Inverse normalization factor */
    } HanningTransformConstants;

    /**
     * @struct RocoFEstimationStates
     * @brief State variables for ROCOF (Rate of Change of Frequency) estimation
     * 
     * Maintains the state needed for ROCOF estimation, including previous frequency
     * values, thresholds, filter coefficients, and filter state.
     */
    typedef struct
    {
        float_p freq_old[NUM_CHANLS];       /**< Previous frequency values for each channel */
        float_p thresholds[3];              /**< ROCOF detection thresholds */
        float_p low_pass_coeff[3];          /**< Low-pass filter coefficients */
        float_p delay_line[NUM_CHANLS][2];  /**< Filter delay line for each channel */
        bool_p state[NUM_CHANLS];           /**< Filter state for each channel */
    } RocoFEstimationStates;

    /**
     * @struct pmu_context
     * @brief PMU estimator context structure
     * 
     * This is the main context structure that encapsulates all state and parameters
     * for a single PMU estimator instance. Multiple independent instances can be
     * created and managed simultaneously.
     * 
     * @note In version 1.7.0 and later, multiple independent instances are supported.
     */
    typedef struct pmu_context
    {
        SynchrophasorEstimatorParams synch_params;  /**< Synchrophasor estimation parameters */
        InternalBuffers buff_params;                /**< Internal buffer pointers */
        HanningTransformConstants hann_params;      /**< Hanning transform constants */
        RocoFEstimationStates rocof_params;         /**< ROCOF estimation state */
        bool_p pmu_initialized;                     /**< Initialization flag */
    } pmu_context;

    /**
     * @brief Initialize a PMU estimator instance
     * 
     * Initializes a PMU estimator context with the specified configuration.
     * The configuration can be loaded from an INI file or passed as a structure.
     * 
     * @param[in,out] ctx Pointer to the PMU context structure to initialize
     * @param[in] cfg Configuration: either a filename (char*) if config_from_ini=1,
     *                or a pointer to an estimator_config structure if config_from_ini=0
     * @param[in] config_from_ini Configuration source flag:
     *                            - CONFIG_FROM_INI (1): Load from INI file
     *                            - CONFIG_FROM_STRUCT (0): Use structure
     * 
     * @return 0 on success, -1 on error
     * 
     * @note The context should not be already initialized. Call pmu_deinit() first
     *       if reinitializing an existing context.
     * 
     * @see pmu_deinit()
     * @see estimator_config
     */
    int pmu_init(pmu_context *ctx, void *cfg, bool_p config_from_ini);

    /**
     * @brief Estimate synchrophasor, frequency, and ROCOF from input signal
     * 
     * This is the main estimation function that processes a window of input samples
     * and produces synchrophasor and ROCOF estimates using the Iterative Enhanced
     * Interpolated DFT algorithm.
     * 
     * @param[in] ctx Pointer to initialized PMU context
     * @param[in] in_signal_windows Pointer to input signal samples. For multi-channel:
     *                              array of pointers to signal windows for each channel.
     *                              Array size must match NUM_CHANLS configuration.
     * @param[in] mid_fracsec Fractional second timestamp (relative to PPS) of the
     *                        midpoint of the observation window. Used for phase correction.
     * @param[out] out_frame Pointer to output PMU frame structure where results are stored
     * 
     * @return 0 on success, -1 on error
     * 
     * @warning The size of in_signal_windows must match the NUM_CHANLS setting used
     *          during compilation. Mismatch will cause segmentation fault.
     * 
     * @note The context must be initialized with pmu_init() before calling this function.
     * 
     * @see pmu_init()
     * @see pmu_frame
     */
    int pmu_estimate(pmu_context *ctx, float_p *in_signal_windows, float_p mid_fracsec, pmu_frame *out_frame);

    /**
     * @brief Deinitialize and free resources of PMU estimator instance
     * 
     * Releases all dynamically allocated memory associated with the PMU context
     * and resets the initialization flag. After calling this function, the context
     * can be reinitialized with pmu_init() if needed.
     * 
     * @param[in,out] ctx Pointer to the PMU context to deinitialize
     * 
     * @return 0 on success, -1 on error
     * 
     * @note Always call this function when finished with a PMU context to prevent
     *       memory leaks.
     * 
     * @see pmu_init()
     */
    int pmu_deinit(pmu_context *ctx);

    /**
     * @brief Print PMU frame data to output stream
     * 
     * Outputs the contents of a pmu_frame structure in a human-readable format
     * to the specified output stream (e.g., stdout, stderr, or a file).
     * 
     * @param[in] frame Pointer to the PMU frame to dump
     * @param[in] stream Output stream (e.g., stdout, stderr, or FILE pointer from fopen)
     * 
     * @return 0 on success, -1 on error
     * 
     * @see pmu_frame
     */
    int pmu_dump_frame(pmu_frame *frame, FILE *stream);

#ifdef __cplusplus
}
#endif

#endif /* PMU_ESTIMATOR_H */