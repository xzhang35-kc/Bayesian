/////////////////////////////////////
// Name: PX-GS-main.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: Parameter Expanded Algorithm

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>

#include "Parameters.h"
#include "DataInput.h"
#include "MCMC.h"
#include "Random.h"
#include "Matrix.h"
#include "MH.h"

gsl_rng * r;

void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, const int & m, const int & N);
void Gibbsampler(gsl_matrix * Y, gsl_matrix * X, gsl_matrix * YY, gsl_matrix ** XX, gsl_matrix * Psigma, gsl_matrix ** Gama, gsl_matrix * JJ,
                 const int & N, const int & P, const int & rep, const int & S, const int & m,
                 gsl_matrix * b, gsl_matrix * C, char * outBetaFile, char * outZFile, char * outRFile, char * outSigmaFile, char * outGamaFile);

int main(int argc, char * argv[])
{
	// define the parameters used in the program
	int N = -1; // number of samples
	int P = -1; // number of covariates
	int rep = -1; // number of repeated measures
	int S = -1; // number of iterations
	int m = -1; // d.f. for the proposed density of correlation
	int m0 = 10; // d.f. for the prior of correlation
	int PD = 1; // prior for Sigma, default is 1 and identity
    int SMH = 1000; // number of iteration for inner MH
    int DFMH = 60; // proposed d.f. for inner MH

	// input/output files
	char * dirName = new char [100];
	strcpy(dirName, "");
	char * xFile = new char [200];
	strcpy(xFile, "");
	char * yFile = new char [200];
	strcpy(yFile, "");
	char * outFilePrefix = new char [200];
	strcpy(outFilePrefix, "");

    // set up parameters
    setPara(argc, argv, N, P, rep, S, m, m0, PD, SMH, DFMH, dirName, yFile, xFile, outFilePrefix);
    printPara(N, P, rep, S, m, m0, PD, SMH, DFMH, dirName, yFile, xFile, outFilePrefix);

    // prepare input data
    gsl_matrix * X = gsl_matrix_alloc(N * rep, P); // covariate matrix N*K by P
    gsl_matrix ** XX = (gsl_matrix **) calloc(N, sizeof(gsl_matrix * ));
	for(int n = 0; n < N; n++)
    {
		XX[n] = gsl_matrix_alloc(rep, P);
    }
    gsl_matrix * Y = gsl_matrix_alloc(N * rep, 1); // response N*K
    gsl_matrix * YY = gsl_matrix_alloc(N, rep);
    DataInput(dirName, yFile, xFile, Y, YY, X, XX, N, rep, P);

    // prepare data
    gsl_matrix * XB = gsl_matrix_alloc(N * rep, 1); // product of X and beta
    gsl_matrix * Psigma = gsl_matrix_alloc(rep, rep);
    gsl_matrix * b = gsl_matrix_alloc(P, 1); // prior mean for beta
    gsl_matrix * C = gsl_matrix_calloc(P, P); // prior covariance for beta

    // find number of cut-off points
    int ordcut;
    FindNumCut(ordcut, Y);
 	gsl_matrix * JJ = gsl_matrix_alloc(rep, 1);
    for(int k = 0; k < rep; k++)
    {
        gsl_matrix_set(JJ, k, 0, ordcut);
    }
    gsl_matrix ** Gama = (gsl_matrix **) calloc(rep, sizeof(gsl_matrix *));
	for(int k = 0; k < rep; k++)
	{
		Gama[k] = gsl_matrix_alloc(gsl_matrix_get(JJ, k, 0), 1);
	}

	// set sigma matrix for the prior distribution of R and D
    setStrSigma(Psigma, b, C, m, m0, PD, 1);

	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

 	// prepare output files
 	char * outDir = new char [200];
 	if(PD == 1) strcpy(outDir, "res-pid-px-gs\\");
    else if(PD == 2) strcpy(outDir, "res-pcs-px-gs\\");
    else  strcpy(outDir, "");
    strcat(dirName, outDir);
    mkdir(dirName);

	char * outBetaFile = new char [200];
	strcpy(outBetaFile, dirName);
	strcat(outBetaFile, "CorrOrdBeta.dat");

    char * outZFile = new char [200];
	strcpy(outZFile, dirName);
	strcat(outZFile, "CorrOrdZ.dat");

    char * outRFile = new char [200];
	strcpy(outRFile, dirName);
	strcat(outRFile, "CorrOrdR.dat");

    char * outSigmaFile = new char [200];
	strcpy(outSigmaFile, dirName);
	strcat(outSigmaFile, "CorrOrdSig.dat");

    char * outGamaFile = new char [200];
	strcpy(outGamaFile, dirName);
    strcat(outGamaFile, "CorrOrdGama.dat");

	// perform Gibbs Sampling
    Gibbsampler(Y, X, YY, XX, Psigma, Gama, JJ, N, P, rep, S, m, b, C, outBetaFile, outZFile, outRFile, outSigmaFile, outGamaFile);

	// free memory;
	gsl_matrix_free(X);
	gsl_matrix_free(Y);
	gsl_matrix_free(YY);
	gsl_matrix_free(XB);
	gsl_matrix_free(Psigma);
	gsl_matrix_free(b);
	gsl_matrix_free(C);
	gsl_matrix_free(JJ);
	for(int n = 0; n < N; n++)
	{
		gsl_matrix_free(XX[n]);
	}
	free(XX);
	for(int k = 0; k < rep; k++)
	{
		gsl_matrix_free(Gama[k]);
	}
	free(Gama);
	delete []dirName;
	delete []yFile;
	delete []xFile;
    delete []outDir;
	delete []outBetaFile;
	delete []outZFile;
	delete []outRFile;
	delete []outSigmaFile;
	delete []outGamaFile;
}

// sample covariance and correlation matrix
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, const int & m, const int & N)
{
	int K = Sigma -> size1;

    gsl_matrix ** WXB, ** T;
    WXB = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
    T = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
    gsl_matrix * Txx = gsl_matrix_calloc(K, K);
    for(int i = 0; i < K; i++)
    {
        for(int j = 0; j < K; j++)
            gsl_matrix_set(Txx, i, j, 0);
    }
	for(int n = 0; n < N; n++)
	{
        WXB[n] = gsl_matrix_alloc(K, 1);
        T[n] = gsl_matrix_alloc(K, K);
    }
    double a, b, c;
    for(int n = 0; n < N; n++)
	{
		for(int k = 0; k < K; k++)
		{
			a = gsl_matrix_get(W, n * K + k, 0);
			//gsl_matrix_set(WW[i], j, 0, a);
			b = gsl_matrix_get(XB, n * K + k, 0);
			//gsl_matrix_set(XXB[i], j, 0, b);
			c = a - b;
			gsl_matrix_set(WXB[n], k, 0, c);
		}

		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, WXB[n], WXB[n], 0.0, T[n]);

		gsl_matrix_add(Txx, T[n]);
    }
    gsl_matrix_add(Txx, Psigma);

    int dfinv = N + m + K + 1;
    InvWishart(Sigma, dfinv, Txx);
    CorrelationM(R, Sigma);

    gsl_matrix_free(Txx);
	for(int n = 0; n < N; n++)
	{
		gsl_matrix_free(WXB[n]);
		gsl_matrix_free(T[n]);
	}
	free(WXB);
	free(T);
}

// Gibbs sampler, main function
void Gibbsampler(gsl_matrix * Y, gsl_matrix * X, gsl_matrix * YY, gsl_matrix ** XX, gsl_matrix * Psigma, gsl_matrix ** Gama, gsl_matrix * JJ,
                 const int & N, const int & P, const int & rep, const int & S, const int & m,
                 gsl_matrix * b, gsl_matrix * C, char * outBetaFile, char * outZFile, char * outRFile, char * outSigmaFile, char * outGamaFile)
{
	gsl_matrix * mu = gsl_matrix_alloc(P, 1);
	gsl_matrix * cov = gsl_matrix_alloc(P, P);
	gsl_matrix * XB = gsl_matrix_alloc(N * rep, 1);

	// here we define Z, beta and R to replace the previous parameters;
	gsl_matrix *newZ = gsl_matrix_alloc(N * rep, 1);
	gsl_matrix *newSig = gsl_matrix_alloc(rep, rep);
	gsl_matrix *newBeta = gsl_matrix_alloc(P, 1);
	gsl_matrix *newR = gsl_matrix_alloc(rep, rep);
	gsl_matrix * SD = gsl_matrix_alloc(rep, rep);

	// set the initial value for R;
	//gsl_matrix_set_identity(RR);
	gsl_matrix_set_identity(newSig);
	gsl_matrix_set_identity(newR);
	gsl_matrix_set_identity(SD);

    gsl_matrix **newGama;
	//initialize newGama
	newGama = (gsl_matrix **) calloc(rep, sizeof(gsl_matrix*));
	for(int k = 0; k < rep; k++)
	{
		newGama[k] = gsl_matrix_calloc(gsl_matrix_get(JJ, k, 0), 1);
		if(gsl_matrix_get(JJ, k, 0) >1 )
			gsl_matrix_set(newGama[k], gsl_matrix_get(JJ, k, 0) -1, 0, 10000.0);
	}
	// set the initial value for Z;
    double a;
    for(int n = 0; n < N * rep; n++)
	{
		a = gsl_matrix_get(Y, n, 0);
		if(a > 0)
			gsl_matrix_set(newZ, n, 0, rltruncnorm(0.0, 1.0, 0.0));
		//else if(fabs(a-999)<0.01)
		//	gsl_matrix_set(newZ, n, 0, gsl_ran_gaussian(r, 1.0));
		else
			gsl_matrix_set(newZ, n, 0, rrtruncnorm(0.0, 1.0, 0.0));
	}

	// open the file for output;
	FILE *betafp, *zfp, *sigfp, *rfp, *gamafp;
	if((betafp = fopen(outBetaFile, "w")) == NULL)
	{
		printf("can not open file %s in Gibbssampler()!\n", outBetaFile);
		exit(1);
	}
	if((zfp = fopen(outZFile, "w")) == NULL)
	{
		printf("can not open file %s in Gibbssampler()!\n", outZFile);
		exit(1);
	}
	if((sigfp = fopen(outSigmaFile,"w")) == NULL)
	{
		printf("can not open file %s in Gibbssampler()!\n", outSigmaFile);
		exit(1);
	}
	if((rfp = fopen(outRFile,"w")) == NULL)
	{
		printf("can not open file %s in Gibbssampler()!\n", outRFile);
		exit(1);
	}
	if((gamafp = fopen(outGamaFile,"w")) == NULL)
	{
		printf("can not open file %s in Gibbssampler()!\n", outGamaFile);
		exit(1);
	}
	// beginning of the loop, compare with the previous routine to know the difference;
	// all matrix indexed by i-1 are replaced by initMatrix;
	// all matrix indexed by i are replaced by newMatrix;
	for(int s = 1; s <= S; s++)
	{
       if((s % 100) == 0)
       {
	     printf("step in Gibbs Sampling : %d\n", s);
       }
       cholbeta(mu, cov, newZ, X, XX, newSig, N, b, C);
       multinorm(newBeta, mu, cov);
       gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, newBeta, 0.0, XB);
       SampleSig(newSig, newR, newZ, XB, Psigma, m, N);
       SampleZ(newZ, newSig, XB, Y, N, newGama, JJ);
       SampleGama(N, rep, 0, YY, newZ, newGama, JJ);
       // this part is for the output;
       gsl_matrix_fprintf(betafp, newBeta, "%g");
       gsl_matrix_fprintf(sigfp, newSig, "%g");
       gsl_matrix_fprintf(rfp, newR, "%g");
       //gsl_matrix_fprintf(gamafp, newGama, "%g");
       for(int k = 0; k < rep; k++)
       {
           gsl_matrix_fprintf(gamafp, newGama[k], "%g");
       }
		//if((i % 100) == 0) gsl_matrix_fprintf(zfp, newZ, "%g");
    }
	// close the file;
	fclose(betafp);
	fclose(zfp);
	fclose(sigfp);
	fclose(rfp);
	fclose(gamafp);
	// free memory;
	gsl_matrix_free(mu);
	gsl_matrix_free(cov);
	gsl_matrix_free(XB);
	gsl_matrix_free(newZ);
	gsl_matrix_free(newR);
	gsl_matrix_free(newBeta);
	gsl_matrix_free(newSig);
	gsl_matrix_free(SD);
	for(int k = 0; k < rep; k++)
	{
		gsl_matrix_free(newGama[k]);
	}
	free(newGama);
}

