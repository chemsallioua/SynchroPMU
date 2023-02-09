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

void dft_r(double* in_ptr, double complex* out_ptr , unsigned int out_len, int n_bins);
double hann(double* out_ptr, unsigned int out_len);
double complex wf(int k, double f, double ampl, double phse, double df, int N,double norm_factor);
double complex D(double k, double N);
double complex whDFT(double k, int N);
int ipDFT(double complex* Xdft, int n_bins, double df, double* amp, double* ph, double* freq);
void e_ipDFT(double complex* Xdft, int n_bins,int window_len, double df, int P, double norm_factor, double* amp, double* ph, double* freq);
void pureTone(double complex* Xpure, int n_bins, double f, double ampl, double phse, double df, int N,double norm_factor);
void print_bins(complex *bins, int n_bins, double df, char* str);
void find3LargestIndx(double arr[], int size, int *km, int *kl, int *kr);

int main() {

    double amp=2;
    double ph=1;
    double freq= 55;
    double ki = 0;
    double fi = 25;
    double n =2048 ;
    double fs = 25600;
    
    double dt = 1/fs;
    double df = fs/n;
    
    int P = 3;
    int Q = 3;
    double signal_window[(int)n];
    double hann_window[(int)n];

    double norm_factor = hann(hann_window, n);

    int i, j;
    for(i=0; i<n; i++){
        signal_window[i] = (amp*cos(2*M_PI*freq*dt*i + ph) + amp*ki*cos(2*M_PI*fi*dt*i + ph))*hann_window[i];
    }

    printf("\n== M-Class Parameters ========================================================\n");
    printf("Signal Fundamental Component | Amp(V): %0.2lf | Ph(rad): %0.2lf | Freq(Hz): %0.2lf\n", amp, ph, freq);
    printf("Interference | I-Mag(%%): %0.2lf | I-Freq(Hz): %0.2lf\n", ki*100, fi);
    printf("------------------------------------------------------------------------------\n");
    printf("Window | SamplingFreq(kS/s): %0.3lf | NCycles: %1.0f | FreqResolution: %0.2lf\n", (float)fs/1000, (n*50/fs), df);
    printf("Iterations | P: %d | Q: %d \n", P, Q);
    printf("===============================================================================\n");


    //Ipdft starts here////////////////////////////////////////////////////////////

    unsigned int n_bins = 11;
    double complex dftbins[n_bins];

    dft_r(signal_window, dftbins, (unsigned int)n , n_bins);

    debug_bins(dftbins, n_bins, df, "Input Signal DFT BINS");

    double freq_f, amp_f, ph_f;
    double freq_i, amp_i, ph_i;
    double complex Xf[n_bins];
    double complex Xi[n_bins];
    double complex Xi_pure[n_bins];

    if(!ipDFT(dftbins, n_bins, df, &amp_f, &ph_f, &freq_f)){
        e_ipDFT(dftbins, n_bins, n, df, P, norm_factor, &amp_f, &ph_f, &freq_f);
    }
        //debug("dft------------------------\n");
    debug("\n[iter-e-ipDFT] ###############################################\n");
    for (i = 0; i < Q; i++)
    {   
        debug("\n[iter-e-ipDFT ITERATION: %d] ------------\n", i+1);

        pureTone(Xf, n_bins, freq_f, amp_f, ph_f, df, n , norm_factor);
        for ( j = 0; j < n_bins; j++){
            Xi[j] = dftbins[j] - Xf[j];
        }
        if(!ipDFT(Xi, n_bins, df, &amp_i, &ph_i, &freq_i)){
            e_ipDFT(Xi, n_bins, n, df, P, norm_factor, &amp_i, &ph_i, &freq_i);
        }
        pureTone(Xi_pure, n_bins, freq_i, amp_i, ph_i, df, n , norm_factor);
        for ( j = 0; j < n_bins; j++)
        {
            Xf[j] = dftbins[j] - Xi_pure[j];
        }
        if(!ipDFT(Xf, n_bins, df, &amp_f, &ph_f, &freq_f)){
            e_ipDFT(Xf, n_bins, n, df, P, norm_factor, &amp_f, &ph_f, &freq_f);
        }
        debug("\nEND iter-e-ipDFT ITERATION --------------------------\n");
    }
    debug("\n[END iter-e-ipDFT] ##############################################\n\n");
    
    //Ipdft finishes here////////////////////////////////////////////////////////////
    
    printf("\n---- [Results] ----------------------------------------------------------\n");
    printf("|\tFREQ: %.10lf | AMP: %.10lf | PH: %.10lf \t|\n", freq_f, 2*amp_f/norm_factor, ph_f);
    printf("-------------------------------------------------------------------------\n");

    return 0;
}

void dft_r(double* in_ptr, double complex* out_ptr , unsigned int input_len, int n_bins){
    // debug("dft------------------------\n");
    int k,n;
    double E = 0;
    double temp_abs;
    for (k = 0 ; k < n_bins ; ++k)
    {
        out_ptr[k] = 0;
        for (n=0 ; n<input_len ; ++n) out_ptr[k] += (in_ptr[n] * cexp(-I*((n * k * M_PI*2 / (double)input_len))));
         
        temp_abs = cabs(out_ptr[k]);
        E += temp_abs*temp_abs;   
    }
}

double hann(double* out_ptr, unsigned int out_len){
    
    double norm_fact =0;
    int i=0;
    for (i=0; i < out_len; i++){
 	   out_ptr[i] = 0.5*(1-cos(2*M_PI*i/out_len));
       norm_fact += out_ptr[i]; 
    }
    return norm_fact;
    
}
int ipDFT(double complex* Xdft, int n_bins, double df, double* amp, double* ph, double* freq){

    int j, k1, k2,k3;
    double Xdft_mag[n_bins]; //magnitude of dft

    debug("\n[ipDFT] ===============================================\n");
    
    debug_bins(Xdft, n_bins, df, "DFT BINS IN ipDFT");

    for(j = 0; j < n_bins; j++){        
        Xdft_mag[j] = cabs(Xdft[j]); 
    }

 
    find3LargestIndx(Xdft_mag, n_bins, &k1, &k2, &k3);

    debug("[%s] k1: %d, k2: %d, k3: %d\n",__FUNCTION__, k1,k2,k3);

    double delta_corr = 2*(Xdft_mag[k3]-Xdft_mag[k2])/(Xdft_mag[k2]+Xdft_mag[k3]+2*Xdft_mag[k1]);

    debug("[%s] delta_corr: %lf\n",__FUNCTION__,delta_corr);
    *freq = (k1+delta_corr)*df;

    if(fabs(delta_corr) <= pow(10,-12)){

        *amp =  Xdft_mag[k1];  
        *ph = carg(Xdft[k1]);
        

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, *freq, *amp, *ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 1; 
    }
    else{
        //ipdft estimated quantities
        *amp = Xdft_mag[k1]*fabs((delta_corr*delta_corr-1)*(M_PI*delta_corr)/sin(M_PI*delta_corr)); 
        *ph = carg(Xdft[k1])-M_PI*delta_corr;
        

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, *freq, *amp, *ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 0;
    }
}
void e_ipDFT(double complex* Xdft, int n_bins,int window_len, double df, int P, double norm_factor, double* amp, double* ph, double* freq){

    //computing the magnitude of the DFT and extracting the largest magnitude and its relative index
    int i,j, p;
    double Amp = *amp;  
    double Phse = *ph;
    double Freq = *freq;
        
    double complex X_neg;
    double complex X_pos[n_bins];
    double X_pos_mag[n_bins];
    int k1, k2,k3, sigma;
    double delta_corr;  

    debug("\n[e_ipDFT] ===============================================\n");

    for(p=0 ; p<P ; p++){ //e-ipdft iterations-------------------------------------
    
        debug("\n[e_ipDFT ITERATION: %d] ------------\n", p+1);

        for(j = 0; j < n_bins; j++){ 
            X_neg = wf(j,-Freq, Amp ,-Phse, df, window_len, norm_factor);           
            X_pos[j] = Xdft[j] - X_neg; 
            X_pos_mag[j] = cabs(X_pos[j]);
        }
        debug_bins(X_pos, n_bins, df, "DFT BINS IN e_ipDFT");

        find3LargestIndx(X_pos_mag, n_bins, &k1, &k2, &k3);
        debug("[%s] k1: %d, k2: %d, k3: %d\n",__FUNCTION__, k1,k2,k3);

        double delta_corr = 2*(X_pos_mag[k3]-X_pos_mag[k2])/(X_pos_mag[k2]+X_pos_mag[k3]+2*X_pos_mag[k1]);
        debug("[%s] delta_corr: %lf\n",__FUNCTION__,delta_corr);
        Freq = (k1+delta_corr)*df;

        if(fabs(delta_corr) <= pow(10,-12)){
        Amp = X_pos_mag[k1];  
        Phse = carg(X_pos[k1]);
        }
        else{
        //estimated quantities
        Amp = X_pos_mag[k1]*fabs((delta_corr*delta_corr-1)*(M_PI*delta_corr)/sin(M_PI*delta_corr)); 
        Phse = carg(X_pos[k1])-M_PI*delta_corr;
        }

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, Freq, Amp, Phse);
        debug("\nEND e_ipDFT ITERATION --------------------------\n");
        
    }
    debug("\n[END e_ipDFT]========================================================\n\n");

    *amp =  Amp;  
    *ph = Phse;
    *freq = Freq;

}
double complex whDFT(double k, int N){
    return -0.25*D(k-1,N) + 0.5*D(k,N) - 0.25*D(k+1,N); 
}
double complex D(double k, double N){
    //debug("k: %f, N:%f, D: %f \n",k, N, cabs(cexp(-I*M_PI*k*(N-1)/N)*sin(M_PI*k)/sin(M_PI*k/N)));
    return cexp(-I*M_PI*k*(N-1)/N)*sin(M_PI*k)/sin(M_PI*k/N);
}
double complex wf(int k, double f, double ampl, double phse, double df, int N,double norm_factor){
    //debug("k: %d, wf: %f \n", k, cabs(ampl*cexp(I*phse)*whDFT(k-(f/df), N)/norm_factor));
    return ampl*cexp(I*phse)*whDFT(k-(f/df), N)/norm_factor;
}
void pureTone(double complex* Xpure, int n_bins, double f, double ampl, double phse, double df, int N,double norm_factor){
    debug("\n[pureTone] ===============================================\n");
    int i;
    for (i = 0; i < n_bins; i++)
    {
      Xpure[i] = wf(i, f, ampl, phse, df, N, norm_factor) + wf(i, -f, ampl, -phse, df, N, norm_factor);
    } 

    
    debug("n_bins: %d | f: %0.3lf | ampl: %0.3lf | phse: %0.3lf | df: %0.3lf | N: %d | norm_factor: %0.3lf\n", n_bins, f, ampl, phse, df, N, norm_factor);
    debug_bins(Xpure, n_bins, df, "DFT BINS PURE TONE");
    debug("[END pureTone] ================================================\n\n");

}

void find3LargestIndx(double arr[], int size, int *km, int *kl, int *kr) {
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

void print_bins(complex *bins, int n_bins, double df, char* str){

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