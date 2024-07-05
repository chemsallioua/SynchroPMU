#ifndef FUNC_STUBS_H
#define FUNC_STUBS_H

/*==============================================================================
  @file func_stubs.h

  Pmu estimator functions stubs.

  Authors: Chemseddine Allioua, Brahim Mazighi

  Copyright (c) 2023.
  All Rights Reserved.
  Confidential and Proprietary - University of Bologna.

  This header file defines macros to stub out arithmetic functions and data types.
  The macros replace the function name with acorresponding "pmue_" function that
  you define, allowing you to implement your own versions of these functions with
  the behavior that you want. You can also stub out a function using an external
  library function, if you include the necessary headers and link with the library.

==============================================================================*/

#ifdef __cplusplus
extern "C"
{
#endif
#include <math.h>
#include <complex.h>

/* Define a macro to stub constants ===========================================*/

/* Define a macro to stub out M_PI */
#define M_PI_p M_PI

/* Define a macro to stub data types ==========================================*/

/* Define a macro to stub out the float data type function */
#define float_p float

/* Define a macro to stub out the complex data type function */
#define complex_p _Complex

/* Define a macro to stub out the complex data type function */
#define uint_p unsigned int

/* Define a macro to stub out the complex data type function */
#define bool_p _Bool

    /* Define a macro to stub function ==========================================*/

    /* Define a macro to stub out the fft for real input function */
    /*
    #define pmue_fft_r(in_ptr, out_ptr, out_len, n_bins) my_fft((in_ptr), (out_ptr), (out_len), (n_bins))

    int_p my_fft(float_p* in_ptr, float_p complex_p* out_ptr , uint_p out_len, uint_p n_bins) {
        //Your implementation of exp here
    }
    */

#ifndef pmue_fft_r
#define pmue_fft_r(in_ptr, out_ptr, out_len, n_bins) dft_r(in_ptr, out_ptr, out_len, n_bins)
#endif

    /* Define a macro to stub out the cabs for complex input */
    /*
    #define pmue_cabs(x) my_cabs(x)

    float_p my_cabs(float_p complex_p x) {
        //Your implementation of exp here
    }
    */

#ifndef pmue_cabs
#define pmue_cabs(x) cabs(x)
#endif

/* Define a macro to stub out the cabs for real input */
/*
#define pmue_fabs(x) my_fabs(x)

float_p my_fabs(float_p x) {
    //Your implementation of exp here
}
*/
#ifndef pmue_fabs
#define pmue_fabs(x) fabs(x)
#endif

/* Define a macro to stub out the carg for real input */
/*
#define pmue_carg(x) my_carg(x)

float_p my_carg(complex_p x) {
    //Your implementation of exp here
}
*/

/* Define a macro to stub out the carg function */
#ifndef pmue_carg
#define pmue_carg(x) carg(x)
#endif

/* Define a macro to stub out the cos function */
/*
#define pmue_cos(x) my_cos(x)

float_p my_cos(float_p x) {
    //Your implementation of exp here
}
*/

/* Define a macro to stub out the cos function */
#ifndef pmue_cos
#define pmue_cos(x) cos(x)
#endif

/* Define a macro to stub out the sin function */
/*
#define pmue_sin(x) my_sin(x)

float_p my_sin(float_p x) {
    //Your implementation of exp here
}
*/

/* Define a macro to stub out the sin function */
#ifndef pmue_sin
#define pmue_sin(x) sin(x)
#endif

/* Define a macro to stub out the sin function */
/*
#define pmue_cexp(x) my_cexp(x)

float_p my_cexp(complex_p x) {
    //Your implementation of exp here
}
*/

/* Define a macro to stub out the cexp function */
#ifndef pmue_cexp
#define pmue_cexp(x) cexp(x)
#endif

#ifdef __cplusplus
}
#endif

#endif /* FUNC_STUBS_H */
