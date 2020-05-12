#pragma once

#define N 3 // Number of states
#define F 6 // Number of features extracted from the signal

#include "dsp.h"
#include "stdlib.h"

class Word {

	char label[10];

	double A[N][N]; // Transition matrix
	double B[N][50];
	double init[N]; // Initial states

	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wc++11-extensions"
	double miu[N][F]; // Mean vector
	double sigma[N][F][F]={0}; // Covariance
	#pragma GCC diagnostic pop

	double alfa[N][50];
	double beta[N][50];


public:

	Word();
	Word(char label[10]);

	double log_Likelihood(double obs[6][50], int nrObservations);

	void forward(double obs[6][50], int nrObservations); 
	void backward(double obs[6][50], int nrObservations); 

	void stateLikelihood(double obs[6][50], int nrObservations);

	void initialStep(double obs[6][50], int nrObservations);
	void emStep(double obs[6][50], int nrObservations); 

	void train(double obs[6][50], int nrObservations);

	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wc++11-extensions"
	double logLikelihood=0;
	#pragma GCC diagnostic pop
	
};


// Vector cu suma elementelor egala cu 1
void normalize(double v[], double len) {
	double sum=0;
	for (int i=0; i<len; i++) {
		sum+=v[i];
	}
	for(int i=0; i<len; i++) {
		v[i]/=sum;
	}
}

// Matrice stocastica -> suma elementelor de pe fiecare rand egala cu 1
void stochasticize(double m[N][N], double size) {
	for(int i=0;i<size;i++) {
		double sum=0;
		for(int j=0; j<size; j++) {
			sum+=m[i][j];
		}
		for(int j=0; j<size; j++) {
			m[i][j]/=sum;
		}
	}
}


#include "Word.cpp"
