#include "gmm.h"



GMM ::GMM() {
	data=NULL;
	miu=length=sigma2=0;
}

GMM ::GMM(int length) {
	this->length=length;
	data=new double[length];
	miu=sigma2=0;
}

void GMM ::setLength(int length) {
	this->length=length;
	if(data==NULL) {
		data=new double[length];
	}
}

void GMM ::setData(double *v) {
	for(int i=0;i<length;i++) {
		data[i]=v[i];
	}
}

void GMM ::calculateMiu() {
	double s=0;
	for(int i=0;i<length;i++) {
		s+=data[i];
	}
	miu = s/(double)length;
    // std::cout<<miu<<" ";
}
void GMM::calculateSigma2() {
	double s=0;
	for(int i=0;i<length;i++) {
		s+=(data[i]-miu)*(data[i]-miu);
	}
    sigma2= s/(double)length;
    // std::cout<<sigma2<<" ";
}

// Functia de densitate de probabilitate 
double GMM:: pdf(double x) {
	return exp(-0.5*(x-miu)*(x-miu)/sigma2)/sqrt(2.0*M_PI*sigma2);
}


