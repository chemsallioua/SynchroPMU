/*==============================================================================
  @file iter_e_ipdft_imp.c

  Source for the implementation of the pmu estimator based on the Iterative
  Enhanced interpolated DFT Algorithm.

  Copyright (c) 2023.
  All Rights Reserved.
  Confidential and Proprietary - University of Bologna.

==============================================================================*/ 

#include "iter_e_ipdft_imp.h"
#include <stdlib.h>

/*GLOBAL VARIABLES DECLARIATION==================*/

// Synchrophasor Estimation Parameters
static unsigned int g_n_channls;        // nuumber of input channels
static unsigned int g_win_len;          // number of samples in observation window
static double g_fs;                     // sample rate in Sample/s
static unsigned int g_n_bins;           // number of bins to define estimation freq band
static unsigned int g_P;                // number of iterations in e_ipDFT
static unsigned int g_Q;                // number of iterations in iter_i_e_ipDFT
static double g_interf_trig;            // trigger value for interference calculation
static double g_df;                     // frequency resolution
static synchrophasor g_phasor;          // global phasor bin
static double g_norm_factor;            // hann normalization factor 

// Dynamically allocated arrays
static double complex* g_Xf;            // Fundamental frequency spectrum bins arry ptr
static double complex* g_Xi;            // Interference frequency spectrum bins arry ptr
static double complex* g_dftbins;       // DFT bins array ptr
static double* g_hann_window;           // hann coefficients array ptr
static double** g_signal_windows;        // input signal windows pointer

// Pmu estimator inizialization flag
static _Bool g_pmu_initialized = 0;

/*STATIC FUNCTIONS PROTOTYPES====================*/

// performs the DFT of a real sampled signal
static double dft_r(double* in_ptr, double complex* out_ptr , unsigned int out_len, unsigned int n_bins);

// helps inizializing the hann coefficients
static double hann(double* out_ptr, unsigned int out_len);

// phasor and frequency estimation main functions 
static void pureTone(double complex* Xpure, synchrophasor phasor);
static int ipDFT(double complex* Xdft, synchrophasor* phasor);
static void e_ipDFT(double complex* Xdft, synchrophasor* phasor);
static void iter_e_ipDFT(complex* dftbins, complex* Xi, complex* Xf, synchrophasor* f_phsr);

// phasor and frequency estimation helping functions
inline static double complex whDFT(double k, int N); 
inline static double complex D(double k, double N);
inline static double complex wf(int k, double f, double ampl, double phse, double df, int N,double norm_factor);
inline static void find3LargestIndx(double arr[], int size, int *km, int *kl, int *kr);

// prints the bins and their index and frequency
static void print_bins(complex *bins, int n_bins, double df, char* str); 

/*FUNCTIONS IMPLEMENTATION====================*/

//pmu initialization function implementation
int pmu_init(void* cfg){

    debug("[%s] Initializing pmu estimator\n",__FUNCTION__);

    estimator_config* config = (estimator_config*)cfg;

    g_win_len = config->win_len;
    g_fs = config->fs;
    g_n_bins = config->n_bins;
    g_P = config->P;
    g_Q = config->Q;
    g_interf_trig = config->interf_trig;
    g_df = g_fs/(double)g_win_len;
    g_n_channls = config->n_chanls;

    if (NULL == (g_Xf = malloc(g_n_bins*sizeof(double complex))  )){
		printf("[%s] ERROR: g_Xf memory allocation failed\n",__FUNCTION__);
		return -1;}
    if (NULL == (g_Xi = malloc(g_n_bins*sizeof(double complex))  )){
		printf("[%s] ERROR: g_Xi memory allocation failed\n",__FUNCTION__);
		return -1;}
    if (NULL == (g_dftbins = malloc(g_n_bins*sizeof(double complex))  )){
		printf("[%s] ERROR: g_Xi memory allocation failed\n",__FUNCTION__);
		return -1;}
    if (NULL == (g_hann_window = malloc(g_win_len*sizeof(double))  )){
		printf("[%s] ERROR: g_hann_window memory allocation failed\n",__FUNCTION__);
		return -1;}

    if (NULL == (g_signal_windows = (double **)malloc(g_n_channls * sizeof(double *))) ){
		printf("[%s] ERROR: g_signal_windows memory allocation failed\n",__FUNCTION__);
		return -1;}

    int i;
    for (i = 0; i < g_n_channls; i++){
        if (NULL == (g_signal_windows[i] = (double *)malloc(g_win_len * sizeof(double)) )){
            printf("[%s] ERROR: g_signal_windows memory allocation failed\n",__FUNCTION__);
            return -1;}
    }
    

    
    g_norm_factor = hann(g_hann_window, g_win_len);

    g_pmu_initialized = 1;

    debug("[%s] Pmu estimator initialized successfully\n",__FUNCTION__);
    return 0;
    
}

//pmu estimation function implementation
int pmu_estimate(double* in_signal_windows[], synchrophasor* out_phasor){

    debug("[%s] pmu_estimate() started\n", __FUNCTION__);

    if(!g_pmu_initialized){
        printf("[%s] ERROR: pmu estimator not initialized, first initialize with pmu_init function\n",__FUNCTION__);
        return -1;
    }

    int i,j, chnl;
    for (chnl = 0; chnl < g_n_channls; chnl++) {
        for(i=0; i<g_win_len; i++){
            g_signal_windows[chnl][i] = in_signal_windows[chnl][i]*g_hann_window[i];
        }
    }

    debug("[%s] windowing on all channels done successfully\n", __FUNCTION__);

    for (chnl = 0; chnl < g_n_channls; chnl++){

        double E = dft_r(g_signal_windows[chnl], g_dftbins, g_win_len , g_n_bins);

        debug_bins(g_dftbins, g_n_bins, g_df, "Input Signal DFT BINS");
        
        e_ipDFT(g_dftbins, &g_phasor);
        pureTone(g_Xf, g_phasor);

        double E_diff = 0;
        
        for ( j = 0; j < g_n_bins; j++){
            g_Xi[j] = g_dftbins[j] - g_Xf[j];
            E_diff += cabs(g_Xi[j]*g_Xi[j]); 
        }

        debug("Energy of Signal Spectrum = %lf | Energy of Difference= %lf\n",E, E_diff);

        if (E_diff > g_interf_trig*E){
            iter_e_ipDFT(g_dftbins, g_Xi, g_Xf, &g_phasor);
        }

        out_phasor[chnl].freq = g_phasor.freq;
        out_phasor[chnl].amp = 2*g_phasor.amp/g_norm_factor;
        out_phasor[chnl].ph = g_phasor.ph;
    }

    return 0;

}

//pmu deinitialization fuunction implementation
int pmu_deinit(){

    debug("[%s] Deinitializing pmu estimator\n",__FUNCTION__);

    if(!g_pmu_initialized){
        printf("[%s] ERROR: pmu estimator not initialized, first initialize with pmu_init function\n",__FUNCTION__);
        return -1;
    }

	free(g_Xf);
	free(g_Xi);
	free(g_dftbins);
    free(g_hann_window);
    free(g_signal_windows);

    int i;
    for (i = 0; i < g_n_channls; i++){
        free(g_signal_windows[i]);
    }
    free(g_signal_windows);

    debug("[%s] Pmu estimator deinitialized successfully\n",__FUNCTION__);

    g_pmu_initialized = 0;
    return 0;
}

static double dft_r(double* in_ptr, double complex* out_ptr , unsigned int out_len, unsigned int n_bins){
    // debug("dft------------------------\n");
    int k,n;
    double E = 0;
    double temp_abs;
    for (k = 0 ; k < n_bins ; ++k)
    {
        out_ptr[k] = 0;
        for (n=0 ; n<out_len ; ++n) out_ptr[k] += (in_ptr[n] * cexp(-I*((n * k * M_PI*2 / (double)out_len))));
         
        temp_abs = cabs(out_ptr[k]);
        E += temp_abs*temp_abs;   
    }
    return E;
}

static double hann(double* out_ptr, unsigned int out_len){
    
    double norm_fact =0;
    int i=0;
    for (i=0; i < out_len; i++){
 	   out_ptr[i] = 0.5*(1-cos(2*M_PI*i/out_len));
       norm_fact += out_ptr[i]; 
    }
    return norm_fact;
    
}

static void pureTone(double complex* Xpure, synchrophasor phasor){
    debug("\n[pureTone] ===============================================\n");
    int i;
    for (i = 0; i < g_n_bins; i++)
    {
      Xpure[i] = wf(i, phasor.freq, phasor.amp, phasor.ph, g_df, g_win_len, g_norm_factor) + wf(i, -phasor.freq, phasor.amp, -phasor.ph, g_df, g_win_len, g_norm_factor);
    } 

    debug("freq: %0.3lf | ampl: %0.3lf | phse: %0.3lf\n", phasor.freq, phasor.amp, phasor.ph);
    debug_bins(Xpure, g_n_bins, g_df, "DFT BINS PURE TONE");
    debug("[END pureTone] ================================================\n\n");

}

static int ipDFT(double complex* Xdft, synchrophasor* phasor){

    int j, k1, k2,k3;
    double Xdft_mag[g_n_bins]; //magnitude of dft

    debug("\n[ipDFT] ===============================================\n");
    
    debug_bins(Xdft, g_n_bins, g_df, "DFT BINS IN ipDFT");

    for(j = 0; j < g_n_bins; j++){        
        Xdft_mag[j] = cabs(Xdft[j]); 
    }

    find3LargestIndx(Xdft_mag, g_n_bins, &k1, &k2, &k3);

    debug("[%s] k1: %d, k2: %d, k3: %d\n",__FUNCTION__, k1,k2,k3);

    double delta_corr = 2*(Xdft_mag[k3]-Xdft_mag[k2])/(Xdft_mag[k2]+Xdft_mag[k3]+2*Xdft_mag[k1]);

    debug("[%s] delta_corr: %lf\n",__FUNCTION__,delta_corr);
    phasor->freq = (k1+delta_corr)*g_df;

    if(fabs(delta_corr) <= pow(10,-12)){

        phasor->amp =  Xdft_mag[k1];  
        phasor->ph = carg(Xdft[k1]);

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, phasor->freq, phasor->amp, phasor->ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 1; 
    }
    else{
        phasor->amp = Xdft_mag[k1]*fabs((delta_corr*delta_corr-1)*(M_PI*delta_corr)/sin(M_PI*delta_corr)); 
        phasor->ph = carg(Xdft[k1])-M_PI*delta_corr;
    
        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, phasor->freq, phasor->amp, phasor->ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 0;
    }
}

static void e_ipDFT(double complex* Xdft, synchrophasor* phasor){

    debug("\n[e_ipDFT] ===============================================\n");

    synchrophasor phsr = *phasor;  
    
    if(!ipDFT(Xdft, &phsr)){        
        int i,j, p;
     
        double complex X_neg;
        double complex X_pos[g_n_bins];

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

    *phasor = phsr;

}

static void iter_e_ipDFT(complex* dftbins, complex* Xi, complex* Xf, synchrophasor* f_phsr){
            
        synchrophasor f_phasor = *f_phsr;
        synchrophasor i_phasor;
        double complex Xi_pure[g_n_bins];
    
        debug("\n[iter-e-ipDFT] ###############################################\n");
        int i,j;
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

inline static double complex whDFT(double k, int N){
    return -0.25*D(k-1,N) + 0.5*D(k,N) - 0.25*D(k+1,N); 
}

inline static double complex D(double k, double N){
    return cexp(-I*M_PI*k*(N-1)/N)*sin(M_PI*k)/sin(M_PI*k/N);
}

inline static double complex wf(int k, double f, double ampl, double phse, double df, int N,double norm_factor){
    return ampl*cexp(I*phse)*whDFT(k-(f/df), N)/norm_factor;
}

inline static void find3LargestIndx(double arr[], int size, int *km, int *kl, int *kr){

  int max_val = -2147483647.0f;
  int max_indx = -1;
  
  for (int i = 0; i < size; i++) {
    if (arr[i] > max_val) {
      max_val = arr[i];
      max_indx = i;
    }
  }
  
  *km = max_indx;
  *kl = max_indx-1;
  *kr = max_indx+1;
}

static void print_bins(complex *bins, int n_bins, double df, char* str){

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
        debug("|%6.1f", cabs(bins[i]));
    }
    debug("|\n---------------------------  ---  --  -\n\n");
}

