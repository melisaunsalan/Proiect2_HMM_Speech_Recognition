#include <iostream>
#include <vector>
#include <complex>
#include <fstream>

#include "dsp.h"
#include "Signal.h"
#include "Word.h"
#include "gmm.h"

using namespace std;


int main() {
	ifstream in("data.txt");

	unsigned int fileNr; // Nr de fisiere audio 
	in>>fileNr;

	Signal s[fileNr]; 

	for(int i=0; i<fileNr; i++) {
		s[i].setFs(8000);
		char name[10];
		in>>name;
		s[i].setLabel(name);
		unsigned int len;
		in>>len;
		s[i].setLength(len);
		for(int j=0;j<len;j++) {
			double val;
			in>>val;
			s[i].setSample(j, val);
		}
	}
	in.close();
  

    for(int i=0;i<fileNr; i++) {
    	s[i].extractFeatures(6,128,32);
    }
    Word word[fileNr];

    for(int i=0;i<fileNr;i++) {
    	std::cout<<s[i].label<<" ";
    	word[i].train(s[i].extractedFeatures,s[i].frameNumber);
    }
   
   
    double logLikelihoods[5][10]; // 60 de fisiere
    // 50 de fisiere pt antrenare si 10 de test
    for(int i=0;i<5;i++) {
    	for(int j=0;j<12;j++) {
    		logLikelihoods[i][j] = word[i*10+j].logLikelihood;
    	}
    }
    
    double testLogLikelihoods[5][2];
    for(int i=0;i<5;i++) {
    	for(int j=0;j<2;j++) {
    		testLogLikelihoods[i][j] = word[50+i*2+j].logLikelihood;
    	}
    }
    
    
    GMM gmm[5];
    for(int i=0;i<5;i++) {
    	gmm[i].setLength(10);
    	gmm[i].setData(logLikelihoods[i]);
    	gmm[i].calculateMiu();
    	gmm[i].calculateSigma2();
    }

    for(int i=0; i<5;i++) {
        for(int j=0;j<10;j++) {
            double prob[]={0};
            double sum=0;
            for(int k=0; k<5;k++) {
                sum+=gmm[k].pdf(logLikelihoods[i][j]);
            }
            for(int k=0;k<5;k++) {
                prob[k]=gmm[k].pdf(logLikelihoods[i][j]);
                prob[k]/=sum;
                // std::cout<<prob[k]<<" ";
            }
            std::cout<<std::endl;
            double pMax=0;
            int kmax=0;
            for(int k=0;k<5;k++) {
                if(prob[k]>pMax) {
                    pMax=prob[k];
                    kmax=k;
                }
            }
            std::cout<<"Guessed: "<<s[kmax*10].label<<" Real: "<<s[i*10+j].label<<std::endl;
        }
        for(int j=0;j<2;j++) {
            double prob[]={0};
            double sum=0;
            for(int k=0; k<5;k++) {
                sum+=gmm[k].pdf(testLogLikelihoods[i][j]);
            }
            for(int k=0;k<5;k++) {
                prob[k]=gmm[k].pdf(testLogLikelihoods[i][j]);
                prob[k]/=sum;
            }
            std::cout<<std::endl;
            double pMax=0;
            int kmax=0;
            for(int k=0;k<5;k++) {
                if(prob[k]>pMax) {
                    pMax=prob[k];
                    kmax=k;
                }
            }
            std::cout<<"Guessed: "<<s[kmax*10].label<<" Real: "<<s[i*10+j].label<<std::endl;
        }
    }
    


    return 0;
}
