
#ifndef ITER_E_IPDFT_IMP_H
#define ITER_E_IPDFT_IMP_H

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

double dft_r(double* in_ptr, double complex* out_ptr , unsigned int out_len, int n_bins);
double hann(double* out_ptr, unsigned int out_len);
double complex wf(int k, double f, double ampl, double phse, double df, int N,double norm_factor);
double complex D(double k, double N);
double complex whDFT(double k, int N);
int ipDFT(double complex* Xdft, int n_bins, double df, double* amp, double* ph, double* freq);
void e_ipDFT(double complex* Xdft, int n_bins,int window_len, double df, int P, double norm_factor, double* amp, double* ph, double* freq);
void iter_e_ipDFT(complex* dftbins, complex* Xi, complex* Xf,double* amp_f,double* ph_f,double* freq_f, int n_bins, int n, double df, int P, int Q, double norm_factor);

void pureTone(double complex* Xpure, int n_bins, double f, double ampl, double phse, double df, int N,double norm_factor);
void print_bins(complex *bins, int n_bins, double df, char* str);

void find_largest_three_indexes(double arr[], int size, int *k1, int *k2, int *k3);

int estimate_synchrophasor(double * signal_window, int n, int fs, int n_bins, int P, int Q, double epsilon,double* out_freq, double* out_amp, double* out_ph);

#endif /* ITER_E_IPDFT_IMP_H */