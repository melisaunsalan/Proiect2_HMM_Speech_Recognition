#include "Word.h"



Word ::Word() {

}


Word ::Word(char label[10]) {
	strcpy(this->label,label);
}


double Word ::log_Likelihood(double obs[6][50], int nrObservations) {
	stateLikelihood(obs, nrObservations);
	// double logLikelihood = 0;
	forward(obs, nrObservations);
	return logLikelihood;
}


void Word ::forward(double obs[6][50], int nrObservations) {
	logLikelihood = 0;

	// Pasul initial
	for(int i=0;i<N;i++) {
		alfa[i][0] = B[i][0] * init[i];
	}

	// Iteratii 
	for(int t=1; t<nrObservations; t++) {

		double temp[N] = {0};
		for(int i=0; i<N; i++) {
			for(int j=0; j<N; j++) {
				temp[i]+=A[j][i]*alfa[j][t-1];
			} 
		}
		for(int i=0;i<N;i++) {
			alfa[i][t] = B[i][t] * temp[i];
		}
		// std::cout<<"alfa \n";
		for(int i=0; i<N; i++) {
			// std::cout<<alfa[i][t]<<" ";
		}
		// std::cout<<std::endl;

		// Scalare
		double sum=0;
		for(int i=0;i<N;i++) {
			sum+=alfa[i][t];
		}
		for(int i=0;i<N;i++) {
			alfa[i][t]/=sum;
		}
		// std::cout<<sum<<std::endl;

		// Update valoare logLikelihood
		logLikelihood+= log(sum);
		// std::cout<<"logLikelihood: "<<logLikelihood<<std::endl;

	}

}


void Word ::backward(double obs[6][50], int nrObservations) {


	// Pasul final
	for(int i=0;i<N;i++) {
		beta[i][nrObservations-1]=1;
	}

	// Iteratii 
	for(int t=nrObservations-2; t>=0; t--) {
		double temp[N] = {0};

		for(int i=0;i<N;i++) {
			temp[i] = beta[i][t+1]*B[i][t+1];
		}

		for(int i=0;i<N;i++) {
			for(int j=0; j<N; j++) {
				beta[i][t]=A[i][j]*temp[j];
			}
		}
		// for(int i=0;i<N;i++) {
		// 	// std::cout<<beta[i][t]<<" ";
		// }

		// Scalare
		double sum=0;
		for(int i=0;i<N;i++) {
			sum+=beta[i][t];
		}
		// std::cout<<"sum= "<<sum<<std::endl;
		for(int i=0;i<N;i++) {
			beta[i][t]/=sum;
		}

	// // std::cout<<"beta \n";
	// 	for(int i=0;i<N;i++) {
	// 		std::cout<<beta[i][t]<<" ";
	// 	}
	// 	std::cout<<std::endl;

	}
}




void Word ::initialStep(double obs[6][50], int nrObservations) {

	// std::cout<<"initialStep called"<<std::endl;

	// Initializarea matricii B cu zerouri
//	B=new double*[N];
//	for(int i=0; i<N; i++) {
//		B[i] = new double[nrObservations];
//		for(int j=0;j<nrObservations; j++) {
//			B[i][j]=0;
//		}
//	}

//	// Initializarea lui alfa cu valori de 0
//	alfa=new double*[N];
//	for(int i=0; i<N; i++) {
//		alfa[i] = new double[nrObservations];
//		for(int j=0;j<nrObservations; j++) {
//			alfa[i][j]=0;
//		}
//	}

//	// Initializarea lui beta cu valori de 0
//	beta=new double*[N];
//	for(int i=0; i<N; i++) {
//		beta[i] = new double[nrObservations];
//		for(int j=0;j<nrObservations; j++) {
//			beta[i][j]=0;
//		}
//	}


	// Initializarea vectorului cu distributia initiala cu valori aleatoare
	for(int i=0; i<N; i++) {
		init[i] = ((double)(rand()%100))/100;
	}
	normalize(init, N);

	// std::cout<<"Initial distribution: "<<std::endl;
	// for(int i=0;i<N;i++) {
	// 	std::cout<<init[i]<<" ";
	// }
	// std::cout<<std::endl;

	// Initializarea matrciii de tranzitie cu valori aleatoare
	for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
			A[i][j]=((double)(rand()%100))/100;
		}
	}
	stochasticize(A,N); 
	// // std::cout<<"Initial A"<<std::endl;
	// for(int i=0;i<N;i++) {
	// 	for(int j=0;j<N;j++) {
	// 		std::cout<<A[i][j]<<" ";
	// 	}
	// 	std::cout<<std::endl;
	// }
	// std::cout<<std::endl;

	// Initializarea matricii de covarianta
	for(int i=0;i<N;i++) {
		covarianceMatrix(obs,nrObservations,sigma[i]);
	}
	for(int i=0;i<N;i++) {
		for(int j=0;j<6;j++) {
			for(int k=0;k<6;k++) {
				if(j!=k) {
					sigma[i][j][k]=0;
				}
			}
		}
	}
	// std::cout<<"Sigma: "<<std::endl;
	// for(int i=0;i<N;i++) {
	// 	for(int j=0;j<6;j++) {
	// 		for(int k=0;k<6;k++) {
	// 			std::cout<<sigma[i][j][k]<<" ";
	// 		}
	// 		std::cout<<std::endl;
	// 	}
	// 	std::cout<<std::endl;
	// }
	// std::cout<<std::endl;
	// std::cout<<determ(sigma[0],6);

	// Initializarea vectorului mediilor
	int indices[N];
	for(int i=0;i<N;i++) {
		indices[i]=rand()%nrObservations;
		for(int j=0;j<i;j++) {
			while(indices[i]==indices[j]) {
				indices[i]=rand()%23;
			}
		}
	}
	for(int i=0;i<N;i++) {
		for(int j=0;j<F;j++) {
			miu[i][j]=obs[j][indices[i]];
		}
	}
	// std::cout<<"Initial miu:"<<std::endl;
	// for(int i=0;i<N;i++) {
	// 	for(int j=0;j<F;j++) {
	// 		std::cout<<miu[i][j]<<" ";
	// 	}
	// 	std::cout<<std::endl;
	// }
	// std::cout<<std::endl;

}



void Word ::stateLikelihood(double obs[6][50], int nrObservations) {

	// Crearea unei distributii normale multivariate
	for(int i=0;i<N;i++) {
		mvnpdf(nrObservations, obs, miu[i], sigma[i], B[i]);
	}
	// std::cout<<std::endl;
	// std::cout<<"Matricea B"<<std::endl;
	// for(int i=0;i<N;i++) {
	// 	for(int j=0;j<nrObservations;j++) {
	// 		std::cout<<B[i][j]<<" ";
	// 	}
	// 	std::cout<<std::endl;
	// }
	// std::cout<<std::endl;
	// std::cout<<std::endl;
}


void Word ::emStep(double obs[6][50], int nrObservations) {
	stateLikelihood(obs, nrObservations);

	forward(obs, nrObservations);
	backward(obs, nrObservations);

	double xiSum[N][N] = {0};
	double gamma[N][nrObservations] ;

	for(int t=0; t<nrObservations-2; t++) {

		double temp[N];
		for(int i=0;i<N;i++) {
			temp[i]=beta[i][t+1]*B[i][t+1];
		}

		double temp2[N][N];
		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				temp2[i][j]=alfa[i][t]*temp[j];
			}
		}

		double temp3[N][N];
		double sum=0;
		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				temp3[i][j] = A[i][j]*temp2[i][j];
				// std::cout<<temp3[i][j];
				sum+=temp3[i][j];
			}
		}
		// std::cout<<"sum: "<<sum<<std::endl;
		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				temp3[i][j]/=(sum+(sum==0));
			}
		}

		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				xiSum[i][j]+=temp3[i][j];
			}
		}

		sum=0;
		for(int i=0; i<N; i++) {
			gamma[i][t]=alfa[i][t]*beta[i][t];
			sum+=gamma[i][t];
		}
		for(int i=0;i<N;i++) {
			gamma[i][t]/=(sum+(sum==0));
		}


	}
	double sum=0;
	for(int i=0; i<N; i++) {
			gamma[i][nrObservations-1]=alfa[i][nrObservations-1]*beta[i][nrObservations-1];
			sum+=gamma[i][nrObservations-1];
		}
	for(int i=0;i<N;i++) {
			gamma[i][nrObservations-1]/=(sum+(sum==0));
	}


	double expectedInit[N];
	for(int i=0;i<N;i++) {
		expectedInit[i]=gamma[i][0];
	}

	double expectedA[N][N];
	for(int i=0;i<N;i++) {
		for(int j=0;j<N;j++) {
			expectedA[i][j]=xiSum[i][j];
		}
	}
	stochasticize(expectedA,N);


	double expectedMiu[N][F] = {0};
	double expectedSigma[N][F][F]={0};


	double gammaStateSum[N]={0};
	for(int i=0;i<N;i++) {
		for(int j=0;j<nrObservations;j++) {
			gammaStateSum[i]+=gamma[i][j];
		}
		if(gammaStateSum[i]==0) {
			gammaStateSum[i]=1;
		}
	}

	double gammaObservations[F][nrObservations];

	for(int s=0;s<N;s++) {
		for (int i = 0; i < F; i++) {
			for (int j=0; j<nrObservations; j++) {
				gammaObservations[i][j] = obs[i][j]*gamma[s][j];
			}
		}

		// Miu
		double sumGammaRow[F]={0};
		for(int i=0;i<F;i++) {
			for(int j=0;j<nrObservations;j++) {
				sumGammaRow[i]+=gammaObservations[i][j];
			}
		}
		for(int i=0;i<F;i++) {
			expectedMiu[s][i]=sumGammaRow[i]/gammaStateSum[s];
		}

		// Sigma

		double tempp[F][F]={0};
		for(int i=0;i<F;i++) {
			for(int j=0;j<F; j++) {
				for(int k=0;k<nrObservations; k++) {
					tempp[i][j]+=gammaObservations[i][k]*obs[j][k];
				}
				tempp[i][j]/=gammaStateSum[s];
				tempp[i][j]-=(miu[s][i]*miu[s][j]);
			}
		}

		// Simetrizarea matricii
		for(int i=0;i<F;i++) {
			for(int j=0;j<F;j++) {
				if(i>j) {
					tempp[i][j]=tempp[j][i];
				}
			}
		}

		// Sigma
		for(int i=0;i<F;i++) {
			for(int j=0;j<F;j++) {
				expectedSigma[s][i][j]=tempp[i][j];
			}
		}

		
	}

	for(int s=0;s<N;s++) {
		for(int i=0;i<F;i++) {
			for(int j=0;j<F;j++) {
				if(i==j) {
					expectedSigma[s][i][j]+=0.01;
				}
			}
		}
	}




	// Urmatorul pas
	// std::cout<<"init \n";
	for(int s=0;s<N;s++) {
		init[s]=expectedInit[s];
		// std::cout<<init[s]<<" ";
	}

	// std::cout<<"miu \n";
	for(int s=0;s<N;s++) {
		for(int i=0;i<F;i++) {
			miu[s][i]=expectedMiu[s][i];
			// std::cout<<miu[s][i]<<" ";
		}
		// std::cout<<"\n";
	}
	// std::cout<<"A \n";
	for(int i=0;i<N;i++) {
		for(int j=0;j<N;j++) {
			A[i][j]=expectedA[i][j];
			// std::cout<<A[i][j]<<" ";
		}
		// std::cout<<"\n";
	}

	// std::cout<<"sigma \n";
	for(int s=0;s<N;s++) {
		for(int i=0;i<F;i++) {
			for(int j=0;j<F;j++) {
				sigma[s][i][j]=expectedSigma[s][i][j];
				// std::cout<<sigma[s][i][j]<<" ";
			}
			// std::cout<<"\n";
		}
		// std::cout<<"\n";
	}
	// std::cout<<"\n";


}

void Word ::train(double obs[6][50], int nrObservations) {

	initialStep(obs, nrObservations);
	log_Likelihood(obs, nrObservations);

	// emStep(obs, nrObservations);
	std::cout<<"logLikelihood: "<<logLikelihood<<std::endl;


}


