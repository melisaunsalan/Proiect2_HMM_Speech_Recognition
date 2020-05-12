#include <cmath>


#pragma once

// Gaussian Mixture Model
class GMM {

public:
	double miu;
	double sigma2;
	double *data;
	int length;

	GMM();
	GMM(int length);
	void setLength(int length);
	void setData(double *v);
	void calculateMiu();
	void calculateSigma2();

	// Functia de densitate de probabilitate
	double pdf(double x);

};

#include "gmm.cpp"