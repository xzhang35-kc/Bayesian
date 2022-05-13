/////////////////////////////////////
// Name: lib-gsm/GSM.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to the GSM sampling method

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

   double a1 = logDet(W);
   double a2 = logDet(sigma);

   Inv(InvW, W);

   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sigma, InvW, 0.0, A);

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
   gsl_matrix_free(InvW);
   gsl_matrix_free(A);
   //compute the Jacobian
   return(0.5 * (K - 1) * t2 - 0.5 * (m0 + K + 1) * a1 - 0.5 * t1 + 0.5 * m0 * a2);
}

// sample SD , standard deviation, square root of variance
void SampleSD(gsl_matrix * SD, gsl_matrix * R, gsl_matrix * Sigma, gsl_matrix * Psigma, const int & m, const int * vt)
{
   int K = Psigma -> size1;
   // double a, b, alpha, beta, c;
   gsl_matrix * InvR = gsl_matrix_alloc(K, K);
   gsl_matrix * DD = gsl_matrix_alloc(K, K);
   gsl_matrix * T = gsl_matrix_alloc(K, K);
   Inv(InvR, R);
   gsl_matrix_set_identity(DD);
   double a, b, c, alpha, beta;
   for (int i = 0; i < K; i++)
   {
      //if(vt[i] == 0)
      //{
      //   a = gsl_matrix_get(Sigma, i, i);
       //  gsl_matrix_set(DD, i, i, sqrt(a));
      //}
      //else if (vt[i] == 1)
     // {
         a = gsl_matrix_get(InvR, i, i);
         b = gsl_matrix_get(Psigma, i, i);
         alpha = m / 2;
         beta = 2 / (a * b);
         c = gsl_ran_gamma(r, alpha, beta);
         gsl_matrix_set(DD, i, i, sqrt(1 / c));
      //}
      //else
      //{
       //  printf("method for nominal data will be implemented later\n");
       //  exit(1);
      //}
   }
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DD, R, 0.0, T);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DD, 0.0, SD);

   // free memory
   free(InvR);
   free(DD);
   free(T);
}

// MH sampling algorithm.
int MHsampler(gsl_matrix * SD, gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Psigma,
   const int & N, const int & mp, const int & m)
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

   int df = mp - K - 1;
   gsl_matrix_memcpy(Wrr, WR);
   gsl_matrix_scale(Wrr, df);

   InvWishart(Wc, mp, Wrr);
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
   double pro1 = Priornew(mp, Wcc, WR);
   double pro2 = Priornew(mp, Wrr, Wc);

   double rho = p1 + pro1 - p2 - pro2;

   if (u <= exp(rho))
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


// adjust Z for continuous, non-missing Y
// Input:
//   sigma: K * K covariance matrix;
//  Y: observed outcomes;
// vt: data type
void AdjustZ(gsl_matrix * Z, gsl_matrix * sigma, gsl_matrix * Y, const int & N, const int * vt)
{
   int K = sigma -> size1;

   int j, w;
   double sigma22, sigma221, u00, a, b, d;
   for (int i = 0; i < N * K; i++)
   {
      j = i % K;
      w = gsl_matrix_get(Y, i, 0) + 0.000001;
      if (w != 999 & vt[j] == 0)
      {
         a = gsl_matrix_get(Y, i, 0);
         b = gsl_matrix_get(sigma, j, j);
         d = a * sqrt(b);
         gsl_matrix_set(Z, i, 0, d);
      }
   }
}


