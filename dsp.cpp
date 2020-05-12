# include "dsp.h"



// Fereastra Hann
 void HannWindow(double w[],unsigned int n){
	for (int i=0; i<n; i++) {
		w[i] = sin(M_PI*i/n);
	}
}


// Transformata Fourier Rapida
std::vector<std::complex<double> > FFT(std::vector<std::complex<double> > &x) {

	int N=x.size();

	if (N==1) {
		return x;
	}

	int M=N/2;

	std::vector<std::complex<double> > Xeven(M,0);
	std::vector<std::complex<double> > Xodd(M,0);

	for(int i=0;i!=M; i++) {
		Xeven[i]=x[2*i];
		Xodd[i]=x[2*i+1];
	}

	std::vector<std::complex<double> > Feven(M,0);
	Feven = FFT(Xeven);
	std::vector<std::complex<double> > Fodd(M,0);
	Fodd = FFT(Xodd);

	std::vector<std::complex<double> > frequencyBins(N,0);
	for(int k=0; k!=N/2; k++) {
		std::complex<double> complexExponential = std::polar(1.0, -2*M_PI*k/N) * Fodd[k];
		frequencyBins[k] = Feven[k] + complexExponential;
		frequencyBins[k+N/2] = Feven[k] - complexExponential;
	}

	return frequencyBins;
}

// Transformata Fourier Inversa
std::vector<std::complex<double> > IFFT(std::vector<std::complex<double> > &X) {

	std::vector<std::complex<double> > x = FFT(X);
	
	for(int i=0; i<X.size(); i++) {
		x[i]/=X.size();
	}
	return x;
}

// Determinantul unei matrici 
double determinant(double m[6][6], int n) {
	int det=0;
	double submatrix[6][6];
	if(n==2) {
		return ((m[0][0]*m[1][1]) - (m[1][0]*m[0][1]));
	}
	else {
		for(int k=0; k<n; k++) {
			int subi = 0;
			for(int i=1;i<n;i++) {
				int subj =0;
				for(int j=0;j<n;j++) {
					if(j==k)
						continue;
					submatrix[subi][subj] =m[i][j];
					subj++;
				}
				subi++;
			}
			det+=(pow(-1,k)*m[0][k]*determinant(submatrix,n-1));
		}
	}
	return det;
}

double determ(double a[6][6],int n) {
  double det=0,  temp[6][6];
  if(n==1) {
    return a[0][0];
  } else if(n==2) {
    det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
    return det;
  } else {
    for(int p=0;p<n;p++) {
      int h = 0;
     int k = 0;
      for(int i=1;i<n;i++) {
        for( int j=0;j<n;j++) {
          if( j==p) {
            continue;
          }
          temp[h][k] = a[i][j];
          k++;
          if(k==n-1) {
            h++;
            k = 0;
          }
        }
      }
      det=det+a[0][p]*pow(-1,p)*determ(temp,n-1);
    }
    return det;
  }
}

// Inversarea unei matrici
void inverseMatrix(double matrix[][20], int order)
{

  double temp;

  // Matricea extinsa
  for (int i = 0; i < order; i++) {
    for (int j = 0; j < 2 * order; j++) {
      if (j == (i + order))
        matrix[i][j] = 1;
    }
  }

  for (int i = order - 1; i > 0; i--) {
     if (matrix[i - 1][0] < matrix[i][0])
     for (int j = 0; j < 2 * order; j++) {
	     temp = matrix[i][j];
	     matrix[i][j] = matrix[i - 1][j];
	     matrix[i - 1][j] = temp;
     }
  }

  for (int i = 0; i < order; i++) {
    for (int j = 0; j < order; j++) {
      if (j != i) {
        temp = matrix[j][i] / matrix[i][i];
        for (int k = 0; k < 2 * order; k++) {
          matrix[j][k] -= matrix[i][k] * temp;
        }
      }
    }
  }

  for (int i = 0; i < order; i++) {

    temp = matrix[i][i];
    for (int j = 0; j < 2 * order; j++) {

      matrix[i][j] = matrix[i][j] / temp;
    }
  }
  return;
}


void vectorXmatrix(double v[6], double m[6][6], double out[6]) {
	for(int i=0;i<6;i++) {
		double s=0;
		for(int j=0;j<6;j++) {
			s+=v[j]*m[j][i];
		}
		out[i]=s;
	}
}

double vectorXvector(double v[6], double w[6]) {
	double s=0;
	for(int i=0;i<6;i++) {
		s+=v[i]*w[i];
	}
	return s;
}

void vectorMvector(double v[6], double w[6], double out[6]) {
	for(int i=0;i<6;i++) {
		out[i]=v[i]-w[i];
	}
}


// Distributia normala multivariata 
void mvnpdf(int n, double x[6][50], double miu[6], double sigma[6][6], double out[n]) {

	double det = determ(sigma,6);

	double invSigma[6][20];
	for(int i=0;i<6;i++) {
		for(int j=0;j<6;j++) {
			invSigma[i][j]=sigma[i][j];
		}
	}
	inverseMatrix(invSigma,6);
	for(int i=0;i<6;i++) {
		for(int j=6;j<12;j++) {
			sigma[i][j-6]=invSigma[i][j];
		}
	}

	for(int i=0;i<n;i++) {
		double xvect[6];
		for(int j=0;j<6;j++) {
			xvect[j] = x[j][i];
		}
		double temp[6];
		vectorMvector(xvect,miu,temp);
		double temp2[6];
		vectorXmatrix(temp,sigma,temp2);
		double exponent = vectorXvector(temp2,temp);
		out[i] = exp(-0.5*exponent)/sqrt(pow(2.0*M_PI, 6) * det);
	}
}

void covarianceMatrix(double obs[6][50], int nrObservations, double cov[6][6]) {

	// Calculul mediei
	double mean[6]={0};
	for(int i=0;i<6;i++) {
		for(int j=0; j<nrObservations;j++) {
			mean[i]+=obs[i][j];
		}
	}

	for(int i=0;i<6;i++) {
		mean[i]/=(double)nrObservations;
	}

	// Calculul covariantei
	for(int i=0; i<6; i++) {
		for(int j=0; j<6; j++) {
			for(int k=0;k<nrObservations;k++) {
				cov[i][j] += (obs[i][k]-mean[i])*(obs[j][k]-mean[j]);
			}
			cov[i][j]/=((double)nrObservations-1);
		}
	}

}
