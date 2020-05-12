#include "Signal.h"



Signal ::Signal() {
	Fs=0;
	length=0;
	data=NULL;
}

Signal ::Signal(int length) {
	this->length=length;
	this->Fs=0;
	this->data = new double[length];
}

Signal ::Signal(int length, int Fs) {
	this->length=length;
	this->Fs=Fs;
	this->data = new double[length];
}

void Signal ::setLabel(char ch[10]) {
	strcpy(this->label, ch);
}

void Signal ::setLength(int length) {
	this->length=length;
	if(this->data==NULL) {
		this->data = new double[length];
	}
}

void Signal ::setFs(int Fs) {
	this->Fs=Fs;
}

double Signal ::getSample(int i) {
	if(i<0 || i>=this->length)
		return 0;
	return data[i];
}

void Signal ::setSample(int i, double newValue) {
	if(i<0 || i>=this->length)
		return;
	data[i]=newValue;
}

void Signal ::print() {
	std::cout<<"( ";
	for(int i=0;i<length-1; i++) {
		std::cout<<data[i]<<", ";
	}
	std::cout<<data[length-1]<<" )"<<std::endl;
}

void Signal ::slice(int start, int stop, double out[]) {
	for(int i=start; i<=stop; i++) {
		out[i-start]=data[i];
	}
}


void Signal ::extractFeatures(unsigned int nrFeatures, unsigned int frameLength, unsigned int overlap) {

	double **buffer;
	buffer = new double*[frameLength]; 
	// extractedFeatures = new double*[nrFeatures];
	int bufferLength = ceil((double)this->length/(double)(frameLength-overlap));
	
	frameNumber=bufferLength;

	for(int i=0; i <frameLength; i++) {
		// extractedFeatures[i] = new double[bufferLength];
		buffer[i] = new double[bufferLength];
	}


	// Crearea buffer-ului cu lungime frameLength si suprapunere overlap
	for(int i=0; i<frameLength; i++) {
		if(i<overlap) {
			buffer[i][0] = 0;
		}
		else {
			buffer[i][0] = data[i-overlap];
		}
	}
	for(int j=1; j<bufferLength-1; j++) {
		for (int i=0; i<overlap; i++) {
			buffer[i][j] = buffer[frameLength - overlap +i][j-1];
		}
		for (int i=overlap; i<frameLength; i++) {
			buffer[i][j] = data[(frameLength-overlap)*j+(i-overlap)];
		}
	}
	for(int i=0; i<frameLength ; i++) {
		if(i<overlap) {
			buffer[i][bufferLength-1] = buffer[frameLength - overlap +i][bufferLength-2];
		}
		else if((i-overlap)<length%(frameLength-overlap)){
			buffer[i][bufferLength-1] = data[(frameLength-overlap)*(bufferLength-1)+(i-overlap)];
		}
		else  {
			buffer[i][bufferLength-1] = 0;
		}
	}


	// Window function
	double window[frameLength];
	HannWindow(window, frameLength);

	// Vector auxiliar pentru operatii cu fiecare frame in parte
	std::vector<std::complex<double> > x(frameLength,0);

	// k parcurge fiecare frame in parte
	for (unsigned int k=0; k<bufferLength; k++) {


		for(unsigned int i=0; i<frameLength; i++) {
			x[i] = window[i]*buffer[i][k];
		}

		// Transformata Fourier a frame-ului inmultit cu fereastra
		std::vector<std::complex<double> > X = FFT(x);
		// for(int i=0;i<frameLength; i++) {
			// std::cout<<X[i].real()<<" + "<<X[i].imag()<<"i ";
		// }

		// Modulul transformatei Fourier
		double fftAbs[frameLength/2];
		for(int i=0;i<frameLength/2;i++) {
			fftAbs[i] = std::abs(X[i]);
		}

		// localMaxima stocheaza pozitiile componentelor spectrale importante
		int localMaxima[frameLength/2]; 
		int maxCounter=0;

		if(fftAbs[0]>fftAbs[1]) {
			localMaxima[maxCounter++]=0;
		}
		for(int i=1;i<frameLength/2-1; i++) {
			if(fftAbs[i]>fftAbs[i-1] && fftAbs[i]>fftAbs[i+1]) {
				localMaxima[maxCounter++]=i;
			}
		}
		if(fftAbs[frameLength/2-1]>fftAbs[frameLength/2-2]) {
			localMaxima[maxCounter++]=frameLength/2-1;
		}

		// Ordonarea pozitiilor varfurilor in ordinea descrescatoare a valorii lor
		for(int i=0;i<maxCounter-1;i++) {
			for(int j=0; j<maxCounter-i-1;j++) {
				if(fftAbs[localMaxima[j]] < fftAbs[localMaxima[j+1]]) {
					double tmp = localMaxima[j];
					localMaxima[j] = localMaxima[j+1];
					localMaxima[j+1] = tmp;
				}
			}
		}



		// for(int i=0;i<nrFeatures;i++) {
		// 	std::cout<<fftAbs[localMaxima[i]]<<" ";
		// }


		for(int i=0; i<nrFeatures; i++) {
			extractedFeatures[i][k]=((double)(localMaxima[i]))/frameLength*Fs;
			// std::cout<<extractedFeatures[i][k]<<" ";
		}

		// std::cout<<std::endl;

	}



}


