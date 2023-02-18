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
#define error(...) fprintf(stderr,__VA_ARGS__)
#else
#define error(...)
#endif
#if LOGGING_LEVEL >= INFO
#define info(...) fprintf(stdout,__VA_ARGS__)
#else
#define info(...)
#endif
#if LOGGING_LEVEL >= DEBUG
#define debug(...) fprintf(stdout,__VA_ARGS__)
#define debug_bins(...) print_bins(__VA_ARGS__)
#else
#define debug(...)
#define debug_bins(...)
#endif

/*MACROS ==================*/

#define D(k, N) (pmue_cexp(-I*M_PI_p*(k)*((N)-1)/(N))*pmue_sin(M_PI_p*(k))/pmue_sin(M_PI_p*(k)/(N)))
#define whDFT(k, N) (-0.25*D((k-1),(N)) + 0.5*D((k),(N)) - 0.25*D((k+1),(N)))
#define wf(k, f, ampl, phse, df, N, norm_factor) (ampl*pmue_cexp(I*phse)*whDFT((k)-((f)/(df)), N)/norm_factor)
#define wrap_angle(rad) (((float_p)(rad) - 2 * M_PI_p * rint((float_p)(rad) / (2 * M_PI_p)))) // wraps angle (rad) in range [-pi; pi]

/*GLOBAL VARIABLES ==================*/

// Synchrophasor Estimation Parameters
static uint_p g_win_len;          // number of samples in observation window
static uint_p g_n_cycles;         // number of cycles at nominal frequency in observation window
static uint_p g_f0;               // fundamental frequency in Hz
static uint_p g_frame_rate;       // frame rate in Frames/second
static uint_p g_fs;               // sample rate in Sample/s
static uint_p g_n_bins;           // number of bins to define estimation freq band
static uint_p g_P;                // number of iterations in e_ipDFT
static uint_p g_Q;                // number of iterations in iter_i_e_ipDFT
static bool_p g_iter_eipdft_enabled;     // flag to enable/disable iterative e_ipDFT
static float_p g_interf_trig;            // trigger value for interference calculation
static float_p g_df;                     // frequency resolution
static phasor g_phasor;                  // global phasor bin
static float_p g_norm_factor;            // hann normalization factor 

// Dynamically allocated arrays
static float_p complex_p* g_Xf;                    // fundamental frequency spectrum bins arry ptr
static float_p complex_p* g_Xi;                    // interference frequency spectrum bins arry ptr
static float_p complex_p* g_dftbins;               // dft bins array ptr
static float_p* g_hann_window;                     // hann coefficients array ptr
static float_p* g_signal_windows[NUM_CHANLS];      // input signal windows pointer

/*ROCOF estimation variables*/
static float_p g_freq_old[NUM_CHANLS];           // represents f(n-1) and it is initialized to f0 Hz
static float_p g_thresholds[3];       	         // thresholds to trigger change in the g_state S1,S2 for rocof estimation
static float_p g_low_pass_coeff[3];              // low pass filter coefficients a1 , b0, b1 respectively in a digital first-order IIR low-pass filter : y[n]=b0x[n]+b1x[n−1]−a1y[n−1]
static float_p g_delay_line[NUM_CHANLS][2];      // represents pre-filter rocof x(n-1) and after filter rocof y(n-1) and it is initialized to zero Hz/s
static bool_p g_state[NUM_CHANLS];               // "0" zero for static conditions, "1" one for dynamic conditions 

// Pmu estimator inizialization flag
static bool_p g_pmu_initialized = 0;

/*STATIC PROTOTYPES ====================*/

// The DFT implementations of a real sampled signal
static int dft_r(float_p* in_ptr, float_p complex_p* out_ptr , uint_p out_len, uint_p n_bins);
static int fft(float_p* in_ptr, float_p complex_p* out_ptr , uint_p out_len);

// helps inizializing the hann coefficients
static float_p hann(float_p* out_ptr, uint_p out_len);

// phasor and frequency estimation main functions 
static void pureTone(float_p complex_p* Xpure, phasor phasor);
static int ipDFT(float_p complex_p* Xdft, phasor* phasor);
static void e_ipDFT(float_p complex_p* Xdft, phasor* out_phasor);
static void iter_e_ipDFT(float_p complex_p* dftbins,float_p complex_p* Xi,float_p complex_p* Xf, phasor* f_phsr);

// phasor and frequency estimation helping functions
inline static void find3LargestIndx(float_p arr[], int size, uint_p *km,uint_p *kl,uint_p *kr);

// pmu estimator configuration functions
static int config_estimator(void* cfg, bool_p config_from_ini);
static int check_config_validity();
static int config_from_file(char* ini_file_name);

// prints the bins and their index and frequency
static void print_bins(float_p complex_p *bins, int n_bins, float_p df, char* str); 

/*IMPLEMENTATION ====================*/

// pmu initialization function implementation
int pmu_init(void* cfg, bool_p config_from_ini){

    // check if the pmu estimator is already initialized
    if(g_pmu_initialized){
        error("[%s] ERROR: pmu estimator already initialized\n",__FUNCTION__);
        return -1;
    }
    if(config_estimator(cfg, config_from_ini)){
        error("[%s] ERROR: pmu estimator configuration failed\n",__FUNCTION__);
        return -1;
    }

    info("[%s] Initializing pmu estimator\n",__FUNCTION__);

    // allocate memory for the global arrays
    if (NULL == (g_Xf = malloc(g_n_bins*sizeof(float_p complex_p))  )){
		error("[%s] ERROR: g_Xf memory allocation failed\n",__FUNCTION__);
		return -1;}
    if (NULL == (g_Xi = malloc(g_n_bins*sizeof(float_p complex_p))  )){
		error("[%s] ERROR: g_Xi memory allocation failed\n",__FUNCTION__);
		return -1;}
    if (NULL == (g_dftbins = malloc(g_win_len*sizeof(float_p complex_p))  )){
		error("[%s] ERROR: g_Xi memory allocation failed\n",__FUNCTION__);
		return -1;}
    if (NULL == (g_hann_window = malloc(g_win_len*sizeof(float_p))  )){
		error("[%s] ERROR: g_hann_window memory allocation failed\n",__FUNCTION__);
		return -1;}
    uint_p i;
    for (i = 0; i < NUM_CHANLS; i++){
        if (NULL == (g_signal_windows[i] = (float_p *)malloc(g_win_len * sizeof(float_p)) )){
            error("[%s] ERROR: g_signal_windows memory allocation failed\n",__FUNCTION__);
            return -1;}
    }

    // initialize the global arrays
    for(i = 0; i < NUM_CHANLS; i++){
        g_delay_line[i][0] = 0.0;
        g_delay_line[i][1] = 0.0;
        g_freq_old[i] = g_f0;
        g_state[i] = 0;
    }
    
    // initialize the hann coefficients
    g_norm_factor = hann(g_hann_window, g_win_len);

    // pmu estimator is initialized successfully
    g_pmu_initialized = 1;

    info("[%s] Pmu estimator initialized successfully\n",__FUNCTION__);
    return 0;
    
}

// pmu estimation function implementation
int pmu_estimate(float_p *in_signal_windows, float_p mid_fracsec ,pmu_frame* out_frame){

    info("[%s] pmu_estimate() started\n", __FUNCTION__);

    // check if pmu estimator is initialized
    if(!g_pmu_initialized){
        error("[%s] ERROR: pmu estimator not initialized, first initialize with pmu_init function\n",__FUNCTION__);
        return -1;
    }

    // input signal windowing
    uint_p j, chnl;
    for (chnl = 0; chnl < NUM_CHANLS; chnl++) {
        for(j=0; j<g_win_len; j++){
            g_signal_windows[chnl][j] = (*((in_signal_windows+chnl*g_win_len) + j))*g_hann_window[j];
        }
    }

    info("[%s] windowing on all channels done successfully\n", __FUNCTION__);

    // synchrophasor estimation for each channel 
    for (chnl = 0; chnl < NUM_CHANLS; chnl++){

        // compuute DFT of input signal
        pmue_fft_r((float_p *)g_signal_windows[chnl], g_dftbins, g_win_len , g_n_bins);

        debug_bins(g_dftbins, g_n_bins, g_df, "Input Signal DFT BINS");
        
        // perform enhanced interpolated DFT
        e_ipDFT(g_dftbins, &g_phasor);
        pureTone(g_Xf, g_phasor);

        float_p E_diff = 0;
        float_p E = 0;

        // compute energy of both input signal and interference frequencies 
        for ( j = 0; j < g_n_bins; j++){

            g_Xi[j] = g_dftbins[j] - g_Xf[j];
            
            if(g_iter_eipdft_enabled){
                E_diff += pmue_cabs(g_Xi[j]*g_Xi[j]);
                E += pmue_cabs(g_dftbins[j]*g_dftbins[j]); 
            }
        }

        debug("Energy of Signal Spectrum (multiplied by epsilon) = %lf | Energy of Difference= %lf\n",E*g_interf_trig, E_diff);

        // check if interference is present to trigger iterative enhanced interpolated DFT
        if(g_iter_eipdft_enabled){
            if (E_diff > g_interf_trig*E){
                iter_e_ipDFT(g_dftbins, g_Xi, g_Xf, &g_phasor);
            }
        }

        // two state rocof estimation
        float_p rocof = (g_phasor.freq - g_freq_old[chnl])*(float_p)g_frame_rate;
        float_p rocof_der = (rocof - g_delay_line[chnl][0])*(float_p)g_frame_rate;

        // update state
        g_state[chnl] = (!g_state[chnl] && ( pmue_fabs(rocof) > g_thresholds[0] || pmue_fabs(rocof_der) > g_thresholds[1] )) ? 1 : g_state[chnl];
        g_state[chnl] = (g_state[chnl] && pmue_fabs(rocof) < g_thresholds[2]) ? 0 : g_state[chnl];
        
        // apply low pass filter if state is 0
        if(!g_state[chnl]){
                out_frame[chnl].rocof = g_low_pass_coeff[1]*rocof + 
                                            g_low_pass_coeff[2]*g_delay_line[chnl][0] - 
                                                g_low_pass_coeff[0]*g_delay_line[chnl][1];
        }else {out_frame[chnl].rocof = rocof;}

        // update old frequency
        g_freq_old[chnl] = g_phasor.freq;

        //update delay line
        g_delay_line[chnl][0] = rocof;    //x(n-1)
        g_delay_line[chnl][1] = out_frame[chnl].rocof; //y(n-1)

        // populate output frame
        out_frame[chnl].synchrophasor.freq = g_phasor.freq;
        out_frame[chnl].synchrophasor.amp = 2*g_phasor.amp/g_norm_factor;
        out_frame[chnl].synchrophasor.ph = wrap_angle(g_phasor.ph - 2*M_PI_p*g_f0*mid_fracsec);

        info("[%s] Estimated on channel: %d\n", __FUNCTION__, chnl);
        
    }

    return 0;

}

// pmu deinitialization fuunction implementation
int pmu_deinit(){

    info("[%s] Deinitializing pmu estimator\n",__FUNCTION__);

    // check if pmu estimator is initialized
    if(!g_pmu_initialized){
        error("[%s] ERROR: pmu estimator not initialized, first initialize with pmu_init function\n",__FUNCTION__);
        return -1;
    }

    // free all allocated memory
	free(g_Xf);
	free(g_Xi);
	free(g_dftbins);
    free(g_hann_window);

    uint_p i;
    for (i = 0; i < NUM_CHANLS; i++){
        free(g_signal_windows[i]);
    }

    info("[%s] Pmu estimator deinitialized successfully\n",__FUNCTION__);

    // set pmu estimator to uninitialized state
    g_pmu_initialized = 0;
    return 0;
}

// dumps quantities of a pmu_frame struct to the specified output stream
int pmu_dump_frame(pmu_frame *frame, FILE *stream){
    
    if (frame == NULL || stream == NULL) {
        error( "Error: NULL pointer passed as argument\n");
        return -1;
    }

    int written = fprintf(stream, "[Synchrophasor] amplitude: %lf, phase: %lf, frequency: %lf, rocof: %lf\n",\
                frame->synchrophasor.amp, wrap_angle(frame->synchrophasor.ph), frame->synchrophasor.freq, frame->rocof);
    if (written < 0) {
        error( "Error: failed to write to stream\n");
        return -1;
    }
    
    return 0;
}

// pmu estimator configuration functions
static int config_estimator(void* cfg, bool_p config_from_ini){

    info("[%s] Configurating pmu estimator\n",__FUNCTION__);

    // pmu estimator parameters initialization from input config
    if(config_from_ini)
    {
        char* ini_file_name = (char*)cfg;
        if(config_from_file(ini_file_name)){
            error("[%s] ERROR: pmu estimator configuration failed\n",__FUNCTION__);
            return -1;
        }
    }
    else
    {
        estimator_config* config = (estimator_config*)cfg;

        g_fs = config->fs;
        g_f0 = config->f0;
        g_n_cycles = config->n_cycles;
        g_frame_rate = config->frame_rate;
        g_n_bins = config->n_bins;
        g_P = config->P;
        g_Q = config->Q;
        g_interf_trig = config->interf_trig;
        g_iter_eipdft_enabled = config->iter_eipdft;

        g_thresholds[0] = config->rocof_thresh[0];
        g_thresholds[1] = config->rocof_thresh[1];
        g_thresholds[2] = config->rocof_thresh[2];

        g_low_pass_coeff[0] = config->rocof_low_pass_coeffs[0];
        g_low_pass_coeff[1] = config->rocof_low_pass_coeffs[1];
        g_low_pass_coeff[2] = config->rocof_low_pass_coeffs[2];
    }

    g_win_len = g_n_cycles*g_fs/g_f0;
    g_df = (float_p)g_fs/(float_p)g_win_len;

    if(check_config_validity()){
        error("[%s] ERROR: pmu estimator configuration failed, config values not valid\n",__FUNCTION__);
        return -1;
    }

    debug("\n[%s] Configuration: g_win_len: %u, g_n_cycles: %u, g_fs: %u, g_f0: %u,\
            g_frame_rate: %u\n g_n_bins: %u, g_P: %u, g_Q: %u, g_interf_trig: %f\n g_df: %f,\
            g_thresholds: %f, %f, %f\n g_low_pass_coeff: %f, %f, %f\n\n",\
            __FUNCTION__,g_win_len, g_n_cycles, g_fs, g_f0, g_frame_rate, g_n_bins, g_P, g_Q,\
            g_interf_trig, g_df, g_thresholds[0], g_thresholds[1], g_thresholds[2],\
            g_low_pass_coeff[0], g_low_pass_coeff[1], g_low_pass_coeff[2]);
    
    info("[%s] Configuration done successfully \n",__FUNCTION__);

    return 0;
}

static int check_config_validity(){

        if(g_f0 != 50 && g_f0 != 60){
            error("[%s] ERROR: nominal frequency: %u not correctly set, (allowed values 50 or 60)\n",__FUNCTION__,g_f0);
            return -1;
        }
    
        if(g_n_cycles <= 0 || g_n_cycles > 50){
            error("[%s] ERROR: window length: %u is not correctly set, must be non-zero, positive and smaller than 50\n",__FUNCTION__, g_n_cycles);
            return -1;
        }
        if(g_fs <= 0 || fmod(g_fs,g_f0) != 0.0){
            error("[%s] ERROR: sample rate: %u is not correctly set, fs must be non-zero, positive and can be devided by f0: %u\n",__FUNCTION__,g_fs, g_f0);
            return -1;
        }

        if(g_frame_rate <= 0 ){
            error("[%s] ERROR: frame rate: %u is not correctly set, must be non-zero and positive\n",__FUNCTION__, g_frame_rate);
            return -1;
        }
        if(g_n_bins < g_n_cycles+2 || g_n_bins > (g_n_cycles*g_fs/g_f0)/2){
            error("[%s] ERROR: number of dft bins: %u for estimation is not correctly set, must be greater than or equal (number of cycles +2): %d and smaller than or equal half the window length: %d \n",__FUNCTION__, g_n_bins, g_n_cycles+2 , (g_n_cycles*g_fs/g_f0)/2);
            return -1;
        }
        if(g_P <= 0){
            error("[%s] ERROR: ipdft iterations: %u is not correctly set, must be non-zero and positive\n",__FUNCTION__, g_P);
            return -1;
        }
        if(g_interf_trig <= 0 || g_interf_trig > 1){
            error("[%s] ERROR: interference threshold: %lf is not correctly set, must be between ]0,1]\n",__FUNCTION__, g_interf_trig);
            return -1;
        }
        if(NUM_CHANLS < 1 || fmod(NUM_CHANLS,1) != 0){
            error("[%s] ERROR: number of channels: %u is not correctly set, must be an integer and at least 1\n",__FUNCTION__, NUM_CHANLS);
            return -1;
        }
    
        if(g_thresholds[0] <= 0){
            error("[%s] ERROR: rocof threshold 1: %lf is not set, must be non-zero and positive\n",__FUNCTION__, g_thresholds[0]);
            return -1;
        }
        if(g_thresholds[1] <= 0){
            error("[%s] ERROR: rocof threshold 2: %lf is not set, must be non-zero and positive\n",__FUNCTION__, g_thresholds[1]);
            return -1;
        }
        if(g_thresholds[2] <= 0){
            error("[%s] ERROR: rocof threshold 3: %lf is not set, must be non-zero and positive\n",__FUNCTION__, g_thresholds[2]);
            return -1;
        }
    
        return 0;

}

static int config_from_file(char* ini_file_name){

    dictionary * ini ;

    ini = iniparser_load(ini_file_name);
    if (ini==NULL) {
        error( "[%s] Error: cannot parse file: %s\n",__FUNCTION__,ini_file_name);
        return -1 ;}
    
    g_n_cycles = (int)iniparser_getint(ini, "signal:n_cycles", 0);
    g_fs = (int)iniparser_getint(ini, "signal:sample_rate", 0);
    g_f0 = (int)iniparser_getint(ini, "signal:nominal_freq", 0);
    g_frame_rate = (int)iniparser_getint(ini, "synchrophasor:frame_rate", 0);
    g_n_bins = (int)iniparser_getint(ini, "synchrophasor:number_of_dft_bins", 0);
    g_P = (int)iniparser_getint(ini, "synchrophasor:ipdft_iterations", 0);
    g_Q = (int)iniparser_getint(ini, "synchrophasor:iter_e_ipdft_iterations", 0);
    g_interf_trig = (float_p)iniparser_getdouble(ini, "synchrophasor:interference_threshold", 0);
    g_iter_eipdft_enabled = (bool_p)iniparser_getboolean(ini, "synchrophasor:iter_e_ipdft_enable", 0);

    g_thresholds[0] = (float_p)iniparser_getdouble(ini, "rocof:threshold_1", 0);
    g_thresholds[1] = (float_p)iniparser_getdouble(ini, "rocof:threshold_2", 0);
    g_thresholds[2] = (float_p)iniparser_getdouble(ini, "rocof:threshold_3", 0);

    g_low_pass_coeff[0] = (float_p)iniparser_getdouble(ini, "rocof:low_pass_filter_1", 0);
    g_low_pass_coeff[1] = (float_p)iniparser_getdouble(ini, "rocof:low_pass_filter_2", 0);
    g_low_pass_coeff[2] = (float_p)iniparser_getdouble(ini, "rocof:low_pass_filter_3", 0);

    iniparser_freedict(ini);

    return 0;
    
}

// performs the DFT of a real sampled signal
static int dft_r(float_p* in_ptr, float_p complex_p* out_ptr , uint_p out_len, uint_p n_bins){
    
    if(((out_len & (out_len - 1)) == 0))
    {   
        float_p g_fft_in[out_len];
        memcpy(g_fft_in, in_ptr, out_len*sizeof(float_p));
        fft(g_fft_in, out_ptr , out_len);

    }else{

        uint_p k,n;

        for (k = 0 ; k < n_bins ; ++k)
        {
            out_ptr[k] = 0;
            for (n=0 ; n<out_len ; ++n) out_ptr[k] += (in_ptr[n] * pmue_cexp(-I*((n * k * M_PI_p*2 / (float_p)out_len))));   
        }
    }

    return 0;
}
static int fft(float_p* in_ptr, float_p complex_p* out_ptr , uint_p out_len){
        
    uint_p half_len = out_len/2;
    uint_p i,k;

    if (out_len == 1) {
        *out_ptr = in_ptr[0];
    } else {
        float_p even[half_len], odd[half_len];
        float_p complex_p even_out[half_len], odd_out[half_len];

        // Split input into even and odd indexed elements
        for (i = 0; i < half_len; i++) {
            even[i] = in_ptr[2*i];
            odd[i] = in_ptr[2*i + 1];
        }

        // Compute FFT of even and odd indexed elements recursively
        fft(even, even_out, half_len);
        fft(odd, odd_out, half_len);

        // Combine results from even and odd indexed elements
        for (k = 0; k < half_len; k++) {
            float_p complex_p t = pmue_cexp(-I * 2 * M_PI * k / out_len) * odd_out[k];
            out_ptr[k] = even_out[k] + t;
            out_ptr[k + half_len] = even_out[k] - t;
        }
    }
    return 0;
}

// helps inizializing the hann coefficients
static float_p hann(float_p* out_ptr, uint_p out_len){
    
    float_p norm_fact =0;
    uint_p i=0;
    for (i=0; i < out_len; i++){
 	   out_ptr[i] = 0.5*(1-pmue_cos(2*M_PI_p*i/out_len));
       norm_fact += out_ptr[i]; 
    }
    return norm_fact;
    
}

// phasor and frequency estimation main functions
static void pureTone(float_p complex_p* Xpure, phasor phasor){
    debug("\n[pureTone] ===============================================\n");
    uint_p i;
    for (i = 0; i < g_n_bins; i++)
    {
      Xpure[i] = wf(i, phasor.freq, phasor.amp, phasor.ph, g_df, g_win_len, g_norm_factor) + wf(i, -phasor.freq, phasor.amp, -phasor.ph, g_df, g_win_len, g_norm_factor);
    } 

    debug("freq: %0.3lf | ampl: %0.3lf | phse: %0.3lf\n", phasor.freq, phasor.amp, phasor.ph);
    debug_bins(Xpure, g_n_bins, g_df, "DFT BINS PURE TONE");
    debug("[END pureTone] ================================================\n\n");

}

static int ipDFT(float_p complex_p* Xdft, phasor* phasor){

    uint_p j, k1, k2,k3;
    float_p Xdft_mag[g_n_bins]; //magnitude of dft

    debug("\n[ipDFT] ===============================================\n");
    
    debug_bins(Xdft, g_n_bins, g_df, "DFT BINS IN ipDFT");

    for(j = 0; j < g_n_bins; j++){        
        Xdft_mag[j] = pmue_cabs(Xdft[j]); 
    }

    find3LargestIndx(Xdft_mag, g_n_bins, &k1, &k2, &k3);

    debug("[%s] k1: %d, k2: %d, k3: %d\n",__FUNCTION__, k1,k2,k3);

    float_p delta_corr = 2*(Xdft_mag[k3]-Xdft_mag[k2])/(Xdft_mag[k2]+Xdft_mag[k3]+2*Xdft_mag[k1]);

    debug("[%s] delta_corr: %lf\n",__FUNCTION__,delta_corr);
    phasor->freq = (k1+delta_corr)*g_df;

    if(pmue_fabs(delta_corr) <= 10e-12){

        phasor->amp =  Xdft_mag[k1];  
        phasor->ph = pmue_carg(Xdft[k1]);

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, phasor->freq, phasor->amp, phasor->ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 1; 
    }
    else{
        phasor->amp = Xdft_mag[k1]*pmue_fabs((delta_corr*delta_corr-1)*(M_PI_p*delta_corr)/pmue_sin(M_PI_p*delta_corr)); 
        phasor->ph = pmue_carg(Xdft[k1])-M_PI_p*delta_corr;
    
        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, phasor->freq, phasor->amp, phasor->ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 0;
    }
}

static void e_ipDFT(float_p complex_p* Xdft, phasor* out_phasor){

    debug("\n[e_ipDFT] ===============================================\n");

    phasor phsr = *out_phasor;  
    
    if(!ipDFT(Xdft, &phsr)){        
        uint_p j, p;
     
        float_p complex_p X_neg;
        float_p complex_p X_pos[g_n_bins];

        for(p=0 ; p<g_P ; p++){ 
        
            debug("\n[e_ipDFT ITERATION: %d] ------------\n", p+1);

            for(j = 0; j < g_n_bins; j++){ 
                X_neg = wf(j,-phsr.freq, phsr.amp ,-phsr.ph, g_df, g_win_len, g_norm_factor);           
                X_pos[j] = Xdft[j] - X_neg; 
            }

            if(ipDFT(X_pos, &phsr)){
                break;
            }
            debug("\nEND e_ipDFT ITERATION --------------------------\n");   
        }
    }
    debug("\n[END e_ipDFT]========================================================\n\n");

    *out_phasor = phsr;

}

static void iter_e_ipDFT(float_p complex_p* dftbins,float_p complex_p* Xi,float_p complex_p* Xf, phasor* f_phsr){
            
        phasor f_phasor = *f_phsr;
        phasor i_phasor;
        float_p complex_p Xi_pure[g_n_bins];
    
        debug("\n[iter-e-ipDFT] ###############################################\n");
        uint_p i,j;
        for (i = 0; i < g_Q; i++)
        {   
            debug("\n[iter-e-ipDFT ITERATION: %d] ------------\n", i+1);


            e_ipDFT(Xi, &i_phasor);
            pureTone(Xi_pure, i_phasor);
            for ( j = 0; j < g_n_bins; j++)
            {
                Xf[j] = dftbins[j] - Xi_pure[j];
            }
            e_ipDFT(Xf, &f_phasor);

            if(i < g_Q-1){

                pureTone(Xf, f_phasor);

                for ( j = 0; j < g_n_bins; j++){
                    Xi[j] = dftbins[j] - Xf[j];
                }
            }

            debug("\nEND iter-e-ipDFT ITERATION --------------------------\n");
        }
        debug("\n[END iter-e-ipDFT] ##############################################\n\n");

        *f_phsr = f_phasor;
}

inline static void find3LargestIndx(float_p arr[], int size, uint_p *km,uint_p *kl,uint_p *kr){

  int max_val = -2147483647.0f;
  int max_indx = -1;
  
  for (int i = 0; i < size; i++) {
    if (arr[i] > max_val) {
      max_val = arr[i];
      max_indx = i;
    }
  }

  if(max_indx == 0 || max_indx == size-1){
    *km = max_indx;
    *kl = max_indx;
    *kr = max_indx;
    return;
  }
  
  *km = max_indx;
  *kl = max_indx-1;
  *kr = max_indx+1;
}

// prints the bins and their index and frequency
static void print_bins(float_p complex_p *bins, int n_bins, float_p df, char* str){

    debug("\n--%s---------------  ---  --  -\n Indx", str);
    for(int i = 0; i < n_bins; i++){
        debug("|%6d", i);
    }
    debug("|\n Freq");
    for(int i = 0; i < n_bins; i++){
        debug("|%6.1f", i*df);
    }
    debug("|\n Bins");
    for(int i = 0; i < n_bins; i++){
        debug("|%6.1f", pmue_cabs(bins[i]));
    }
    debug("|\n---------------------------  ---  --  -\n\n");
}
