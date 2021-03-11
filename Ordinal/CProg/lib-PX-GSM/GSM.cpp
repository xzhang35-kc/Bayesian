/////////////////////////////////////
// Name: GSM.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: routines for GSM method

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>

#include "matrix.h"
#include "Random.h"

extern gsl_rng * r;


double Priornew(const int & m0, gsl_matrix * sigma, gsl_matrix * W)
{
	int K = W -> size1;

	gsl_matrix * InvW = gsl_matrix_alloc(K, K);
    gsl_matrix * A = gsl_matrix_alloc(K, K);

	double a1 = Deter(W);
	double a2 = Deter(sigma);

	Inv(InvW, W);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sigma, InvW, 0.0, A);

	double t1 = 0;
    double t2 = 0;
    double b1, b2;
	for(int i = 0; i < K; i++)
	{
		b1 = gsl_matrix_get(A, i, i);
		b2 = gsl_matrix_get(W, i, i);
		t1 = t1 + b1;
		t2 = t2 + log(b2);
    }

	//free memory
	gsl_matrix_free(InvW);
	gsl_matrix_free(A);
	//compute the Jacobian
	return(0.5 * (K - 1) * t2 - 0.5 * (m0 + K + 1) * a1 - 0.5 * t1 + 0.5 * m0 * a2);
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

// MH sampling algorithm.
int MHsampler(gsl_matrix * SD, gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
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


