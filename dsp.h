# pragma once 

#include <cmath>
#include <vector>
#include <complex>
#include "dsp.cpp"
#include "Signal.h"



// Fereastra Hann
void HannWindow(double w[], unsigned int n);

// Transformata Fourier Rapida
std::vector<std::complex<double> > FFT(std::vector<std::complex<double> > &x);
// Transformata Fourier Inversa
std::vector<std::complex<double> > IFFT(std::vector<std::complex<double> > &X);

// Determinantul unei matrici
double determinant(double m[6][6], int n);

// Inversa unei matrici 
void inverseMatrix(double matrix[][20], int order);

// Distributia normala multivariata
void mvnpdf(int n, double x[6][50], double miu[6], double sigma[6][6], double out[n]);

// Matricea de covarianta
void covarianceMatrix(double obs[6][50], int nrObservations, double cov[6][6]);

