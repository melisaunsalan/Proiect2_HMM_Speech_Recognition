#pragma once

#include <cmath>
#include <iostream>
#include <cstring>
#include "dsp.h"



class Signal {
	
	int Fs;

public:
	// double **extractedFeatures;
	double *data;
	int length;
	double extractedFeatures[6][50];
	int frameNumber;
	char label[10];

public:
	Signal();
	Signal(int length);
	Signal(int length, int Fs);
	void setLabel(char ch[10]);
	void setLength(int length);
	void setFs(int Fs);
	double getSample(int i);
	void setSample(int i, double newValue);
	void print();
	void slice(int start, int stop, double out[]);
	void extractFeatures(unsigned int F, unsigned int frameLength, unsigned int overlap);

};

#include "Signal.cpp"

