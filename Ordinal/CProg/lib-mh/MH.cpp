/////////////////////////////////////
// Name: lib-mh/MH.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to the MH method

#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>

#include "Matrix.h"
#include "Random.h"

extern gsl_rng * r;


// obtain covariance matrix RR from WR for mixed data which used in the posterior density calculation
/*
void Transform(gsl_matrix * RR, gsl_matrix * WR, gsl_matrix * JJ)
{
   int K = WR -> size1;
   gsl_matrix * R = gsl_matrix_alloc(K, K);
   gsl_matrix * DI = gsl_matrix_alloc(K, K);
   gsl_matrix * DI1 = gsl_matrix_alloc(K, K);
   gsl_matrix_set_identity(DI);
   CorrelationM(R, WR);
   double a;
   int b;
   for (int i = 0; i < K; i++)
   {
      b = gsl_matrix_get(JJ, i, 0) + 0.0000001;
      if (b <= 1)
      {
         a = gsl_matrix_get(WR, i, i);
         gsl_matrix_set(DI, i, i, sqrt(a));
      }
   }

   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, DI, R, 0.0, DI1);
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, DI1, DI, 0.0, RR);
   symmetric(RR);

   //free of memory
   gsl_matrix_free(DI);
   gsl_matrix_free(DI1);
   gsl_matrix_free(R);
}
*/

// obtain covariance matrix RR from WR for mixed data which used in the posterior density calculation
void Transform(gsl_matrix * RR, gsl_matrix * WR, const int * vt)
{
	double a;
	int K = WR -> size1;
	gsl_matrix * R = gsl_matrix_alloc(K, K);
	gsl_matrix * DI = gsl_matrix_alloc(K, K);
	gsl_matrix * DI1 = gsl_matrix_alloc(K, K);
	gsl_matrix_set_identity(DI);

	for(int i = 0; i < K; i++)
	{
		if(vt[i] == 1)
		{
			a = gsl_matrix_get(WR, i, i);
			gsl_matrix_set(DI, i, i, 1.0/sqrt(a));
		}
		else if(vt[i] == 2)
      {
         printf("method for nominal data will be implemented later\n");
      }
	}

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DI,  WR, 0.0, DI1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DI1, DI, 0.0, RR);

	//free of memory
	gsl_matrix_free(DI);
	gsl_matrix_free(DI1);
	gsl_matrix_free(R);
}


// compute the prior for both R and D,
// in fact, in this subroutine, we calculate the density of Wishart time the Jacobian
// without the consideration of the constant.
double Priornew(const int & m0, gsl_matrix * sigma, gsl_matrix * W)
{
   int K = W -> size1;

   gsl_matrix * Invsigma = gsl_matrix_alloc(K, K);
   gsl_matrix * A = gsl_matrix_alloc(K, K);

   double a1 = 0.5 * logDet(W);
   double a2 = 0.5 * logDet(sigma);

   Inv(Invsigma, sigma);

   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Invsigma, W, 0.0, A);

   double t1 = 0;
   double t2 = 0;
   double b1, b2;
   for (int i = 0; i < K; i++)
   {
      b1 = gsl_matrix_get(A, i, i);
      b2 = gsl_matrix_get(W, i, i);
      t1 = t1 + b1;
      t2 = t2 + log(b2);
   }

   //free memory
   gsl_matrix_free(Invsigma);
   gsl_matrix_free(A);

   return(0.5 * (K - 1) * t2 + (m0 - K - 1) * a1 - 0.5 * t1 - m0 * a2);
}

// compute log posterior density of correlation matrix not including the prior part.
double Posterior(const int & N, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB)
{
   int K = R -> size1;
   gsl_matrix * Q1 = gsl_matrix_alloc (N * K, 1);
   gsl_matrix * Q2 = gsl_matrix_alloc (K, K);
   gsl_matrix * Sig = gsl_matrix_calloc (K, K);
   gsl_matrix * Siginv = gsl_matrix_alloc (K, K);
   gsl_matrix ** QQ;
   QQ = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
   for (int i = 0; i < N; i++)
   {
      QQ[i] = gsl_matrix_alloc(K, 1);
   }

   gsl_matrix * T1 = gsl_matrix_alloc(1, K);
   gsl_matrix * T2 = gsl_matrix_alloc(1, 1);

   gsl_matrix_memcpy(Q1, Z);
   gsl_matrix_sub (Q1, XB);

   gsl_matrix_memcpy(Q2, R);

   // choleskey decomposition
   gsl_linalg_cholesky_decomp(Q2);

   // obtain lower triangle matrix Sig
   for (int i = 0; i < K; i++)
   {
      for (int j = 0; j <= i; j++)
      {
         gsl_matrix_set (Sig, i, j, gsl_matrix_get (Q2, i, j));
      }
   }

   Inv(Siginv, Sig);
   double a;
   double Txx = 0;
   for (int i = 0; i < N; i++)
   {
      for (int j = 0; j < K; j++)
      {
         a = gsl_matrix_get(Q1, i * K + j, 0);
         gsl_matrix_set(QQ[i], j, 0, a);
      }

      gsl_blas_dgemm (CblasTrans, CblasTrans, 1.0, QQ[i], Siginv, 0.0, T1);

      gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, T1, T1, 0.0, T2);

      Txx = Txx + gsl_matrix_get(T2, 0, 0);
   }

   double sum1 = (-N) * 0.5 * logDet(R);
   double sum2 = (-1.0 / 2.0) * Txx;

   // free memory;
   gsl_matrix_free(Q1);
   gsl_matrix_free(Q2);
   gsl_matrix_free(Sig);
   gsl_matrix_free(Siginv);
   gsl_matrix_free(T1);
   gsl_matrix_free(T2);
   for (int i = 0; i < N; i++)
   {
      gsl_matrix_free(QQ[i]);
   }
   free(QQ);
   return(sum1 + sum2);
}

// M-H sampling algorithm.
// The information for input is described in the main function.
int MHsampler(gsl_matrix * WR, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
           const int & N, const int & m, const int & m0, gsl_matrix * JJ, const int * vt)
{
   int K = WR -> size1;
   gsl_matrix * Wc = gsl_matrix_alloc(K, K);
   gsl_matrix * Wcc = gsl_matrix_alloc(K, K);
   gsl_matrix * Wrr = gsl_matrix_alloc(K, K);
   gsl_matrix * RRc = gsl_matrix_alloc(K, K);
   gsl_matrix * RR0 = gsl_matrix_alloc(K, K);

   gsl_matrix_memcpy(Wrr, WR);
   gsl_matrix_scale(Wrr, 1.0 / m);

   Wishart(Wc, m, Wrr);

   gsl_matrix_memcpy(Wcc, Wc);
   gsl_matrix_scale(Wcc, 1.0 / m);

   Transform(RRc, Wc, vt);
   Transform(RR0, WR, vt);

   double u = gsl_rng_uniform (r);

   double p1 = Priornew(m0, Psigma, Wc);
   double p2 = Priornew(m0, Psigma, WR);
   double post1 = Posterior(N, RRc, Z, XB);
   double post2 = Posterior(N, RR0, Z, XB);
   double pro1 = Priornew(m, Wcc, WR);
   double pro2 = Priornew(m, Wrr, Wc);

   double rho = p1 + post1 + pro1 - p2 - post2 - pro2;

   if (u <= exp(rho))
   {
      gsl_matrix_memcpy(WR, Wc);
   }
   gsl_matrix_free(Wc);
   gsl_matrix_free(Wcc);
   gsl_matrix_free(Wrr);
   gsl_matrix_free(RR0);
   gsl_matrix_free(RRc);
   if (u <= exp(rho))
   {
      return(1);
   }
   else
   {
      return(0);
   }
}

