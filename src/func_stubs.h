/**
 * @file func_stubs.h
 * @brief Function and type stubs for PMU estimator
 * 
 * This header file provides a flexible stubbing mechanism for arithmetic functions
 * and data types used in the PMU estimator. It allows users to replace standard
 * library implementations with custom versions or external library functions.
 * 
 * The stubbing system enables:
 * - Custom implementations of mathematical functions
 * - Integration with hardware-accelerated libraries
 * - Type substitution for different platforms or precision requirements
 * - Easy testing and mocking of mathematical operations
 * 
 * @section usage Usage
 * 
 * To use a custom implementation, define the macro before including this header:
 * @code
 * #define pmue_sin(x) my_custom_sin(x)
 * #include "func_stubs.h"
 * @endcode
 * 
 * Or define your function and use it directly:
 * @code
 * float_p my_custom_fft(float_p* in_ptr, float_p complex_p* out_ptr, 
 *                       uint_p out_len, uint_p n_bins) {
 *     // Your custom FFT implementation
 * }
 * #define pmue_fft_r(in_ptr, out_ptr, out_len, n_bins) \
 *         my_custom_fft((in_ptr), (out_ptr), (out_len), (n_bins))
 * @endcode
 * 
 * @author Chemseddine Allioua, Brahim Mazighi
 * @copyright Copyright (c) 2023. All Rights Reserved.
 *            Confidential and Proprietary - University of Bologna.
 */

#ifndef FUNC_STUBS_H
#define FUNC_STUBS_H

#ifdef __cplusplus
extern "C"
{
#endif
#include <math.h>
#include <complex.h>

/* Define a macro to stub constants ===========================================*/

/** @brief Mathematical constant PI (can be overridden) */
#define M_PI_p M_PI

/* Define a macro to stub data types ==========================================*/

/** @brief Floating-point type stub (default: float) */
#define float_p float

/** @brief Complex number type stub (default: _Complex) */
#define complex_p _Complex

/** @brief Unsigned integer type stub (default: unsigned int) */
#define uint_p unsigned int

/** @brief Boolean type stub (default: _Bool) */
#define bool_p _Bool

    /* Define a macro to stub function ==========================================*/

    /**
     * @brief FFT for real input stub
     * @param in_ptr Input array pointer (real values)
     * @param out_ptr Output array pointer (complex values)
     * @param out_len Output array length
     * @param n_bins Number of frequency bins
     * @return Status code
     * 
     * Example custom implementation:
     * @code
     * #define pmue_fft_r(in_ptr, out_ptr, out_len, n_bins) \
     *         my_fft((in_ptr), (out_ptr), (out_len), (n_bins))
     * 
     * int my_fft(float_p* in_ptr, float_p complex_p* out_ptr, 
     *            uint_p out_len, uint_p n_bins) {
     *     // Your FFT implementation here
     * }
     * @endcode
     */
#ifndef pmue_fft_r
#define pmue_fft_r(in_ptr, out_ptr, out_len, n_bins) dft_r(in_ptr, out_ptr, out_len, n_bins)
#endif

    /**
     * @brief Complex absolute value (magnitude) stub
     * @param x Complex number
     * @return Magnitude of complex number
     * 
     * Can be replaced with custom implementation:
     * @code
     * #define pmue_cabs(x) my_cabs(x)
     * 
     * float_p my_cabs(float_p complex_p x) {
     *     // Your implementation here
     * }
     * @endcode
     */
#ifndef pmue_cabs
#define pmue_cabs(x) cabs(x)
#endif

    /**
     * @brief Floating-point absolute value stub
     * @param x Real number
     * @return Absolute value
     * 
     * Can be replaced with custom implementation:
     * @code
     * #define pmue_fabs(x) my_fabs(x)
     * 
     * float_p my_fabs(float_p x) {
     *     // Your implementation here
     * }
     * @endcode
     */
#ifndef pmue_fabs
#define pmue_fabs(x) fabs(x)
#endif

    /**
     * @brief Complex argument (phase angle) stub
     * @param x Complex number
     * @return Phase angle in radians
     * 
     * Can be replaced with custom implementation:
     * @code
     * #define pmue_carg(x) my_carg(x)
     * 
     * float_p my_carg(complex_p x) {
     *     // Your implementation here
     * }
     * @endcode
     */
#ifndef pmue_carg
#define pmue_carg(x) carg(x)
#endif

    /**
     * @brief Cosine function stub
     * @param x Angle in radians
     * @return Cosine of x
     * 
     * Can be replaced with custom implementation:
     * @code
     * #define pmue_cos(x) my_cos(x)
     * 
     * float_p my_cos(float_p x) {
     *     // Your implementation here
     * }
     * @endcode
     */
#ifndef pmue_cos
#define pmue_cos(x) cos(x)
#endif

    /**
     * @brief Sine function stub
     * @param x Angle in radians
     * @return Sine of x
     * 
     * Can be replaced with custom implementation:
     * @code
     * #define pmue_sin(x) my_sin(x)
     * 
     * float_p my_sin(float_p x) {
     *     // Your implementation here
     * }
     * @endcode
     */
#ifndef pmue_sin
#define pmue_sin(x) sin(x)
#endif

    /**
     * @brief Complex exponential function stub
     * @param x Complex number
     * @return e raised to the power of x
     * 
     * Can be replaced with custom implementation:
     * @code
     * #define pmue_cexp(x) my_cexp(x)
     * 
     * float_p complex_p my_cexp(complex_p x) {
     *     // Your implementation here
     * }
     * @endcode
     */
#ifndef pmue_cexp
#define pmue_cexp(x) cexp(x)
#endif

#ifdef __cplusplus
}
#endif

#endif /* FUNC_STUBS_H */
