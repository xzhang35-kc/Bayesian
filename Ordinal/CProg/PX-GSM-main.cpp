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

void SampleSD(gsl_matrix * SD, gsl_matrix * R, gsl_matrix * Psigma, const int & m);
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, const int & m, const int & N);
void Gibbsampler(gsl_matrix * Y, gsl_matrix * X, gsl_matrix * YY, gsl_matrix ** XX, gsl_matrix * Psigma, gsl_matrix ** Gama, gsl_matrix * JJ,
                 const int & N, const int & P, const int & rep, const int & S, const int & m, const int & PD, const int & SMH, const int &DFMH,
                 gsl_matrix * b, gsl_matrix * C, char * outBetaFile, char * outZFile, char * outRFile, char * outSigmaFile, char * outGamaFile);
int MHsamplerGSM(gsl_matrix * SD, gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
              const int & N, const int & pm, const int & m);

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
    setStrSigma(Psigma, b, C, m, m0, PD, 2);

	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

 	// prepare output files
 	char * outDir = new char [200];
 	if(PD == 1) strcpy(outDir, "res-pid-px-gsm\\");
    else if(PD == 2) strcpy(outDir, "res-pcs-px-gsm\\");
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
    Gibbsampler(Y, X, YY, XX, Psigma, Gama, JJ, N, P, rep, S, m, PD, SMH, DFMH,
                b, C, outBetaFile, outZFile, outRFile, outSigmaFile, outGamaFile);

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

// sample SD , standard deviation, square root of variance
void SampleSD(gsl_matrix * SD, gsl_matrix * R, gsl_matrix * Psigma, const int & m)
{
	int K = Psigma -> size1;
//	double a, b, alpha, beta, c;

    gsl_matrix * InvR = gsl_matrix_alloc(K, K);
    gsl_matrix * DD = gsl_matrix_alloc(K, K);
    gsl_matrix * T = gsl_matrix_alloc(K, K);
    Inv(InvR, R);
    gsl_matrix_set_identity(DD);
    double a, b, c, alpha, beta;
    for(int i = 0; i < K; i++)
    {
        a = gsl_matrix_get(InvR, i, i);
        b = gsl_matrix_get(Psigma, i, i);
        alpha = m / 2;
        beta = 2 / (a * b);
        c = gsl_ran_gamma(r, alpha, beta);
        gsl_matrix_set(DD, i, i, sqrt(1 / c));
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DD, R, 0.0, T);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DD, 0.0, SD);

    // free memory
	free(InvR);
	free(DD);
	free(T);
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
                 const int & N, const int & P, const int & rep, const int & S, const int & m, const int & PD, const int & SMH, const int &DFMH,
                 gsl_matrix * b, gsl_matrix * C, char * outBetaFile, char * outZFile, char * outRFile,  char * outSigmaFile, char * outGamaFile)
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
    gsl_matrix * WR = gsl_matrix_alloc(rep, rep);

	// set the initial value for R;
	gsl_matrix_set_identity(newSig);
	gsl_matrix_set_identity(newR);
	if(PD == 1)
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
	int sp = 0;
	for(int s = 1; s <= S; s++)
	{
	    if((s % 100) == 0)
        {
            printf("step in Gibbs Sampling : %d\n", s);
            if(PD == 2)
                printf("number accepted in mh: %d\n", sp);
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, newBeta, 0.0, XB);

        if(PD == 1)
        {
            SampleSD(SD, newR, Psigma, m);
        }
       else if(PD == 2)
       {
           gsl_matrix_memcpy(SD, newSig);
           gsl_matrix_memcpy(WR, newSig);
           sp = 0;
           for(int i = 0; i < SMH; i++)
                sp = sp + MHsamplerGSM(SD, WR, newR, newZ, XB, Psigma, N, DFMH, m);
       }

       SampleZ(newZ, SD, XB, Y, N, newGama, JJ);
       SampleGama(N, rep, 0, YY, newZ, newGama, JJ);
       SampleSig(newSig, newR, newZ, XB, Psigma, m, N);
       cholbeta(mu, cov, newZ, X, XX, newSig, N, b, C);
       multinorm(newBeta, mu, cov);

       // this part is for the output;
       gsl_matrix_fprintf(betafp, newBeta, "%g");
       gsl_matrix_fprintf(sigfp, newSig, "%g");
       gsl_matrix_fprintf(rfp, newR, "%g");

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
    gsl_matrix_free(WR);
}

// MH sampling algorithm.
int MHsamplerGSM(gsl_matrix * SD, gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
              const int & N, const int & pm, const int & m)
{
	int K = R -> size1;
	gsl_matrix * Wc = gsl_matrix_alloc(K, K);
	gsl_matrix * Wcc = gsl_matrix_alloc(K, K);
	gsl_matrix * Wrr = gsl_matrix_alloc(K, K);
	gsl_matrix * Dc = gsl_matrix_alloc(K, K);
	gsl_matrix * D = gsl_matrix_alloc(K, K);
	gsl_matrix * DcRDc = gsl_matrix_alloc(K, K);
	gsl_matrix * DcR = gsl_matrix_alloc(K, K);
	gsl_matrix * DRD = gsl_matrix_alloc(K, K);
	gsl_matrix * DR = gsl_matrix_alloc(K, K);
	gsl_matrix_set_identity(Dc);
    gsl_matrix_set_identity(D);

	int df = pm - K - 1;
	gsl_matrix_memcpy(Wrr, WR);
	gsl_matrix_scale(Wrr, df);

    InvWishart(Wc, pm, Wrr);
	gsl_matrix_memcpy(Wcc, Wc);
	gsl_matrix_scale(Wcc, df);

	double a, b;
	for(int i = 0; i < K; i++)
	{
	    a = gsl_matrix_get(Wc, i, i);
		gsl_matrix_set(Dc, i, i, sqrt(a));
		b = gsl_matrix_get(WR, i, i);
		gsl_matrix_set(D, i, i, sqrt(b));
    }

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D, R, 0.0, DR);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DR, D, 0.0, DRD);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Dc, R, 0.0, DcR);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DcR, Dc, 0.0, DcRDc);

    double u = gsl_rng_uniform(r);

	double p1 = Priornew(m, Psigma, DcRDc);
	double p2 = Priornew(m, Psigma, DRD);
	double pro1 = Priornew(pm, Wcc, WR);
	double pro2 = Priornew(pm, Wrr, Wc);

    double rho = p1 + pro1 - p2 - pro2;

	if(u <= exp(rho))
	{
		gsl_matrix_memcpy(WR, Wc);
		gsl_matrix_memcpy(SD, DcRDc);
		gsl_matrix_free(Wc);
        gsl_matrix_free(Wcc);
        gsl_matrix_free(Wrr);
        gsl_matrix_free(Dc);
        gsl_matrix_free(D);
        gsl_matrix_free(DR);
        gsl_matrix_free(DRD);
        gsl_matrix_free(DcR);
        gsl_matrix_free(DcRDc);
		return(1);
	}
	gsl_matrix_free(Wc);
	gsl_matrix_free(Wcc);
	gsl_matrix_free(Wrr);
	gsl_matrix_free(Dc);
	gsl_matrix_free(D);
	gsl_matrix_free(DR);
	gsl_matrix_free(DRD);
	gsl_matrix_free(DcR);
	gsl_matrix_free(DcRDc);
	return(0);
}

