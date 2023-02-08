#include "iter_e_ipdft_imp.h"

int estimate_synchrophasor(double * signal_window, int n, int fs, int n_bins, int P, int Q, double epsilon,double* out_freq, double* out_amp, double* out_ph){

    double df = fs/n;

    double hann_window[(int)n];
    double complex dftbins[n_bins];

    double norm_factor = hann(hann_window, n);
    double E = dft_r(signal_window, dftbins, (unsigned int)n , n_bins);

    debug_bins(dftbins, n_bins, df, "Input Signal DFT BINS");

    double freq_f, amp_f, ph_f;
    double complex Xf[n_bins];
    double complex Xi[n_bins];
    
    e_ipDFT(dftbins, n_bins, n, df, P, norm_factor, &amp_f, &ph_f, &freq_f);
    pureTone(Xf, n_bins, freq_f, amp_f, ph_f, df, n , norm_factor);

    double E_diff = 0;
    int j;
    for ( j = 0; j < n_bins; j++){
        Xi[j] = dftbins[j] - Xf[j];
        E_diff += cabs(Xi[j]*Xi[j]); 
    }

    debug("Energy of Signal Spectrum = %lf | Energy of Difference= %lf\n",E, E_diff);

    if (E_diff > epsilon*E){
        iter_e_ipDFT(dftbins, Xi, Xf, &amp_f, &ph_f, &freq_f, n_bins, n, df, P, Q, norm_factor);
    }

    *out_freq = freq_f;
    *out_amp = 2*amp_f/norm_factor;
    *out_ph = ph_f;

    return 0;

}
double dft_r(double* in_ptr, double complex* out_ptr , unsigned int input_len, int n_bins){
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
    return E;
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

    find_largest_three_indexes(Xdft_mag, n_bins, &k1, &k2, &k3);

    debug("[%s] k1: %d, k2: %d, k3: %d\n",__FUNCTION__, k1,k2,k3);

    double delta_corr = 2*(Xdft_mag[k3]-Xdft_mag[k2])/(Xdft_mag[k2]+Xdft_mag[k3]+2*Xdft_mag[k1]);

    debug("[%s] delta_corr: %lf\n",__FUNCTION__,delta_corr);

    if(fabs(delta_corr) <= pow(10,-12)){

        *amp =  Xdft_mag[k1];  
        *ph = carg(Xdft[k1]);
        *freq = k1*df;

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, *freq, *amp, *ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 1; 
    }
    else{
        //ipdft estimated quantities
        *amp = Xdft_mag[k1]*fabs((delta_corr*delta_corr-1)*(M_PI*delta_corr)/sin(M_PI*delta_corr)); 
        *ph = carg(Xdft[k1])-M_PI*delta_corr;
        *freq = (k1+delta_corr)*df;

        debug("[%s] freq: %.10lf, amp (not normalized): %.3lf, ph: %.3lf\n",__FUNCTION__, *freq, *amp, *ph);
        debug("\n[END ipDFT] ===============================================\n\n");

        return 0;
    }
}
void e_ipDFT(double complex* Xdft, int n_bins,int window_len, double df, int P, double norm_factor, double* amp, double* ph, double* freq){

    debug("\n[e_ipDFT] ===============================================\n");

    double Amp = *amp;  
    double Phse = *ph;
    double Freq = *freq;
    
    if(!ipDFT(Xdft, n_bins, df, &Amp, &Phse, &Freq)){        
        //computing the magnitude of the DFT and extracting the largest magnitude and its relative index
        int i,j, p;

            
        double complex X_neg;
        double complex X_pos[n_bins];
        double X_pos_mag[n_bins];
        int k1, k2,k3, sigma;
        double delta_corr;  

        

        for(p=0 ; p<P ; p++){ //e-ipdft iterations-------------------------------------
        
            debug("\n[e_ipDFT ITERATION: %d] ------------\n", p+1);

            for(j = 0; j < n_bins; j++){ 
                X_neg = wf(j,-Freq, Amp ,-Phse, df, window_len, norm_factor);           
                X_pos[j] = Xdft[j] - X_neg; 
                //X_pos_mag[j] = cabs(X_pos[j]);
            }

            if(ipDFT(X_pos, n_bins, df, &Amp, &Phse, &Freq)){
                break;
            }
            debug("\nEND e_ipDFT ITERATION --------------------------\n");   
        }
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

void find_largest_three_indexes(double arr[], int size, int *k1, int *k2, int *k3) {
  int i;
  float first, second, third;
  int first_index, second_index, third_index;
  first = second = third = -2147483647.0f;
  
  for (i = 0; i < size; i++) {
    if (arr[i] > first) {
      third = second;
      second = first;
      first = arr[i];
      third_index = second_index;
      second_index = first_index;
      first_index = i;
    } else if (arr[i] > second) {
      third = second;
      second = arr[i];
      third_index = second_index;
      second_index = i;
    } else if (arr[i] > third) {
      third = arr[i];
      third_index = i;
    }
  }
  
  *k1 = first_index;
  *k2 = second_index;
  *k3 = third_index;
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

void iter_e_ipDFT(complex* dftbins, complex* Xi, complex* Xf,double* amp_f,double* ph_f,double* freq_f, int n_bins, int n, double df, int P, int Q, double norm_factor){
            
        double Amp  =  *amp_f;
        double Phse =  *ph_f;
        double Freq =  *freq_f;

        double freq_i, amp_i, ph_i;
        double complex Xi_pure[n_bins];
    
        debug("\n[iter-e-ipDFT] ###############################################\n");
        int i,j;
        for (i = 0; i < Q; i++)
        {   
            debug("\n[iter-e-ipDFT ITERATION: %d] ------------\n", i+1);


            e_ipDFT(Xi, n_bins, n, df, P, norm_factor, &amp_i, &ph_i, &freq_i);
            pureTone(Xi_pure, n_bins, freq_i, amp_i, ph_i, df, n , norm_factor);
            for ( j = 0; j < n_bins; j++)
            {
                Xf[j] = dftbins[j] - Xi_pure[j];
            }
            e_ipDFT(Xf, n_bins, n, df, P, norm_factor, &Amp, &Phse, &Freq);

            if(i < Q-1){

                pureTone(Xf, n_bins, Freq, Amp, Phse, df, n , norm_factor);

                for ( j = 0; j < n_bins; j++){
                    Xi[j] = dftbins[j] - Xf[j];
                }
            }

            debug("\nEND iter-e-ipDFT ITERATION --------------------------\n");
        }
        debug("\n[END iter-e-ipDFT] ##############################################\n\n");

        *amp_f = Amp;
        *ph_f = Phse;
        *freq_f = Freq;
}