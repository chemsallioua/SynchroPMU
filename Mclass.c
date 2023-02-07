#include <stdio.h>
#include <complex.h>
#include <math.h>

void dft_r(double* in_ptr, double complex* out_ptr , unsigned int out_len, int n_bins);
double hann(double* out_ptr, unsigned int out_len);
double complex wf(int k, double f, double ampl, double phse, double df, int N,double norm_factor);
double complex D(double k, double N);
double complex whDFT(double k, int N);
int ipDFT(double complex* Xdft, int n_bins, double df, double* amp, double* ph, double* freq);
void e_ipDFT(double complex* Xdft, int n_bins,int window_len, double df, int P, double norm_factor, double* amp, double* ph, double* freq);
void pureTone(double complex* Xpure, int n_bins, double f, double ampl, double phse, double df, int N,double norm_factor);

void find_largest_three_indexes(float arr[], int size, int *k1, int *k2, int *k3);

int main() {

    double amp=2;
    double ph=1;
    double freq= 50;
    double ki = 0.1;
    double fi = 25;
    double n =5000 ;
    double fs = n/5*50;
    
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

    //Ipdft starts here////////////////////////////////////////////////////////////

    unsigned int n_bins = 11;
    double complex dftbins[n_bins];

    dft_r(signal_window, dftbins, (unsigned int)n , n_bins);

    double freq_f, amp_f, ph_f;
    double freq_i, amp_i, ph_i;
    double complex Xf[n_bins];
    double complex Xi[n_bins];
    double complex Xi_pure[n_bins];

    if(!ipDFT(dftbins, n_bins, df, &amp_f, &ph_f, &freq_f)){
        e_ipDFT(dftbins, n_bins, n, df, P, norm_factor, &amp_f, &ph_f, &freq_f);
    }
    pureTone(Xf, n_bins, freq_f, amp_f, ph_f, df, n , norm_factor);
    for ( j = 0; j < n_bins; j++){
        Xi[j] = dftbins[j] - Xf[j];
    }
    if(!ipDFT(Xi, n_bins, df, &amp_i, &ph_i, &freq_i)){ 
        //printf("dft------------------------\n");
        for (i = 0; i < Q; i++)
        {   
            e_ipDFT(Xi, n_bins, n, df, P, norm_factor, &amp_i, &ph_i, &freq_i);
            pureTone(Xi_pure, n_bins, freq_i, amp_i, ph_i, df, n , norm_factor);
            for ( j = 0; j < n_bins; j++)
            {
                Xf[j] = dftbins[j] - Xi_pure[j];
            }
            if(!ipDFT(Xf, n_bins, df, &amp_f, &ph_f, &freq_f)){
                e_ipDFT(Xf, n_bins, n, df, P, norm_factor, &amp_f, &ph_f, &freq_f);
            }
        }
    }

    //Ipdft finishes here////////////////////////////////////////////////////////////

    printf("freq: %.10lf, amp: %.10lf, ph: %.10lf\n", freq_f, 2*amp_f/norm_factor, ph_f);

    return 0;
}

void dft_r(double* in_ptr, double complex* out_ptr , unsigned int input_len, int n_bins){
    // printf("dft------------------------\n");
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

    int j, k1, k2;
    double Xdft_mag[n_bins]; //magnitude of dft

    for(j = 0; j < n_bins; j++){        
        Xdft_mag[j] = cabs(Xdft[j]); 
    }

    find2LargstIndx(Xdft_mag, n_bins,&k1,&k2);

    int sigma = ( k2 > k1 ) ? 1:-1; //sign of the delta correction
    double delta_corr = sigma*(2*Xdft_mag[k1+sigma]-Xdft_mag[k1])/(Xdft_mag[k1+sigma]+Xdft_mag[k1]);
    if(delta_corr == 0){

        *amp =  Xdft_mag[k1];  
        *ph = carg(Xdft[k1]);
        *freq = k1*df;

        return 1; 
    }
    else{
        //ipdft estimated quantities
        *amp = Xdft_mag[k1]*fabs((delta_corr*delta_corr-1)*(M_PI*delta_corr)/sin(M_PI*delta_corr)); 
        *ph = carg(Xdft[k1])-M_PI*delta_corr;
        *freq = (k1+delta_corr)*df;

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
    int k1, k2, sigma;
    double delta_corr;  

    for(p=0 ; p<P ; p++){ //e-ipdft iterations-------------------------------------
        for(j = 0; j < n_bins; j++){ 
            X_neg = wf(j,-Freq, Amp ,-Phse, df, window_len, norm_factor);           
            X_pos[j] = Xdft[j] - X_neg; 
            X_pos_mag[j] = cabs(X_pos[j]);
        }
        find2LargstIndx(X_pos_mag, n_bins,&k1,&k2);

        sigma = ( k2 > k1 ) ? 1:-1; 
        delta_corr =sigma*(2*X_pos_mag[k1+sigma]-X_pos_mag[k1])/(X_pos_mag[k1+sigma]+X_pos_mag[k1]);
        Amp = X_pos_mag[k1]*fabs((delta_corr*delta_corr-1)*(M_PI*delta_corr)/sin(M_PI*delta_corr)); 
        Phse = carg(X_pos[k1])-M_PI*delta_corr;
        Freq = (k1+delta_corr)*df;
    }
    *amp =  Amp;  
    *ph = Phse;
    *freq = Freq;

}
double complex whDFT(double k, int N){
    return -0.25*D(k-1,N) + 0.5*D(k,N) - 0.25*D(k+1,N); 
}
double complex D(double k, double N){
    return cexp(-I*M_PI*k*(N-1)/N)*sin(M_PI*k)/sin(M_PI*k/N);
}
double complex wf(int k, double f, double ampl, double phse, double df, int N,double norm_factor){
    return ampl*cexp(I*phse)*whDFT(k-(f/df), N)/norm_factor;
}
void pureTone(double complex* Xpure, int n_bins, double f, double ampl, double phse, double df, int N,double norm_factor){
    int i;
    for (i = 0; i < n_bins; i++)
    {
      Xpure[i] = wf(i, f, ampl, phse, df, N, norm_factor) + wf(i, -f, ampl, -phse, df, N, norm_factor);
    }  
}

    

void find_largest_three_indexes(float arr[], int size, int *k1, int *k2, int *k3) {
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