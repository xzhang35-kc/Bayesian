/////////////////////////////////////
// Name: libsrc/MCMC.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to the MCMC sampling

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

#include "Matrix.h"
#include "Random.h"

extern gsl_rng * r;


// obtain the mean vector of "mu" and co-variance matrix "cov" of beta
// input latent variable Z, covariate X, sigma, and sample size N
void cholbeta(gsl_matrix * mu, gsl_matrix * cov, gsl_matrix * Z, gsl_matrix ** XX,
              gsl_matrix * Sigma, const int & N, gsl_matrix * b, gsl_matrix * InvC)
{
   symmetric(Sigma);
   int K = Sigma -> size1;
   int P = InvC -> size1;
   gsl_matrix ** ZZ;
   ZZ = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
   for (int i = 0; i < N; i++)
   {
      ZZ[i] = gsl_matrix_alloc(K, 1);
   }
   gsl_matrix * Q = gsl_matrix_alloc(K, K);
   gsl_matrix * Sig = gsl_matrix_calloc(K, K);
   gsl_matrix * Siginv = gsl_matrix_alloc(K, K);
   gsl_matrix * T1 = gsl_matrix_alloc(P, K);
   gsl_matrix * T2 = gsl_matrix_alloc(P, P);
   gsl_matrix * T3 = gsl_matrix_alloc(K, 1);
   gsl_matrix * T4 = gsl_matrix_alloc(P, 1);
   gsl_matrix * Txx = gsl_matrix_calloc(P, P);
   gsl_matrix * Txz = gsl_matrix_calloc(P, 1);
   gsl_matrix * CIb = gsl_matrix_calloc(P, 1);

   // copy sigma to Q, so sigma will not be changed;
   gsl_matrix_memcpy(Q, Sigma);

   // choleskey decomposition
   gsl_linalg_cholesky_decomp(Q);

   // get the lower triangle matrix Sig
   double a;
   for (int i = 0; i < K; i++)
   {
      for (int j = 0; j <= i; j++)
      {
         a = gsl_matrix_get(Q, i, j);
         gsl_matrix_set(Sig, i, j, a);
      }
   }
   Inv(Siginv, Sig);

   for (int i = 0; i < N; i++)
   {
      for (int j = 0; j < K; j++)
      {
         a = gsl_matrix_get(Z, i * K + j,0);
         gsl_matrix_set(ZZ[i], j, 0, a);
      }

      gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, XX[i], Siginv, 0.0, T1);

      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T1, T1, 0.0, T2);

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Siginv, ZZ[i], 0.0, T3);

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T1, T3, 0.0, T4);

      gsl_matrix_add(Txx, T2);

      gsl_matrix_add(Txz, T4);
   }

   gsl_matrix_add(Txx, InvC);

   Inv(cov, Txx);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, InvC, b, 0.0, CIb);

    gsl_matrix_add(Txz, CIb);

   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cov, Txz, 0.0, mu);

    // free of memory;
   gsl_matrix_free(Sig);
   gsl_matrix_free(Siginv);
   gsl_matrix_free(Q);
   gsl_matrix_free(T1);
   gsl_matrix_free(T2);
   gsl_matrix_free(T3);
   gsl_matrix_free(T4);
   gsl_matrix_free(Txx);
   gsl_matrix_free(Txz);
   gsl_matrix_free(CIb);
   for (int i = 0; i < N; i++)
   {
      gsl_matrix_free(ZZ[i]);
   }
   free(ZZ);
}

// obtain the sample of latent variable Z
// Input:
//   sigma: K * K covariance matrix;
//   XB: (N * K) * 1 matrix;
//  Y: observed discrete outcomes;
//  N: sample size.
void SampleZ(gsl_matrix * Z, gsl_matrix * sigma, gsl_matrix * XB, gsl_matrix * Y,
             const int & N, gsl_matrix ** Gama, gsl_matrix * JJ, const int * vt)
{
   int K = sigma -> size1;
   gsl_matrix * sigma11 = gsl_matrix_alloc (K - 1, K - 1);
   gsl_matrix * Isigma11 = gsl_matrix_alloc (K - 1, K - 1);
   gsl_matrix * sigma21 = gsl_matrix_alloc (K - 1, 1);
   gsl_matrix * T1 = gsl_matrix_alloc (1, K - 1);
   gsl_matrix * T2 = gsl_matrix_alloc (1, 1);
   gsl_matrix * T3 = gsl_matrix_alloc (K - 1, 1);
   gsl_matrix * T4 = gsl_matrix_alloc (K - 1,1);
   gsl_matrix * T5 = gsl_matrix_alloc (1, 1);

   int j, k, n, w, num;
   double sigma22, sigma221, u00, a, b;
   for (int i = 0; i < N * K; i++)
   {
      j = i % K;
      k = (i - j) / K;
      sigma22 = gsl_matrix_get(sigma, j, j);
      Mat_I_J(sigma11, sigma, j);
      Inv(Isigma11, sigma11);
      Mat_IJ(sigma21, sigma, j);

      gsl_blas_dgemm (CblasTrans,CblasNoTrans, 1.0, sigma21, Isigma11, 0.0, T1);
      gsl_blas_dgemm (CblasNoTrans,CblasNoTrans, 1.0, T1, sigma21, 0.0, T2);
      b = gsl_matrix_get(T2, 0, 0);
      sigma221 = sqrt(sigma22 - b);
      for (n = 0; n < K - 1; n++)
      {
         if ((k * K + n) < i)
         {
            a = gsl_matrix_get(Z, k * K + n, 0);
            gsl_matrix_set(T3, n, 0, a);
         }
         else
         {
            a = gsl_matrix_get(Z, k * K + n + 1, 0);
            gsl_matrix_set(T3, n, 0, a);
         }
      }
      for (n = 0; n < K - 1; n++)
      {
         if ((k * K + n) < i)
         {
            a = gsl_matrix_get(XB, k * K + n, 0);
            gsl_matrix_set(T4, n, 0, a);
         }
         else
         {
            a = gsl_matrix_get(XB, k * K + n + 1, 0);
            gsl_matrix_set(T4, n, 0, a);
         }
      }

      gsl_matrix_sub (T3, T4);

      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0, T1, T3, 0.0, T5);

      u00 = gsl_matrix_get(XB, i, 0) + gsl_matrix_get(T5, 0, 0);

      w = gsl_matrix_get(Y, i, 0) + 0.000001;
      num = gsl_matrix_get(JJ, j, 0) + 0.000001;
      if (w == 999)
      {
         a = gsl_ran_gaussian (r, sigma221) + u00;
         gsl_matrix_set(Z, i, 0, a);
      }
      else if (vt[j] == 1) // for ordinal variable
      {
         if (w == 0)
         {
            a = rrtruncnorm(u00, sigma221, 0.0);
            gsl_matrix_set(Z, i, 0, a);
         }
         else if (w == num - 1)
         {
            a = rltruncnorm(u00, sigma221, gsl_matrix_get(Gama[j], num - 2, 0));
            gsl_matrix_set(Z, i, 0, a);
         }
         else
         {
            for (n = 1; n < num - 1; n++)
            {
               if (w == n)
               {
                  a = truncnorm(u00, sigma221, gsl_matrix_get(Gama[j], n - 1, 0), gsl_matrix_get(Gama[j], n, 0));
                  gsl_matrix_set(Z, i, 0, a);
               }
            }
         }
      }
      else if(vt[j] == 2)
      {
         printf("method for nominal data will be implemented later\n");
         exit(1);
      }
   }

   // free matrix;
   gsl_matrix_free(sigma11);
   gsl_matrix_free(Isigma11);
   gsl_matrix_free(sigma21);
   gsl_matrix_free(T1);
   gsl_matrix_free(T2);
   gsl_matrix_free(T3);
   gsl_matrix_free(T4);
   gsl_matrix_free(T5);
}

// sample the Gama matrix (cut-off point)
void SampleGama(const int & N, const int & K, gsl_matrix * YY, gsl_matrix * Z, gsl_matrix ** Gama, gsl_matrix * JJ, const int * vt)
{
   double a;
   gsl_matrix * ZZ = gsl_matrix_alloc(N, K);
   for (int i = 0; i < N; i++)
   {
      for (int j = 0; j < K; j++)
      {
         a = gsl_matrix_get(Z, i * K + j, 0);
         gsl_matrix_set(ZZ, i, j, a);
      }
   }

   double b, t1, z;
   int ll, yy;
   for (int i = 0; i < K; i++)
   {
      if(vt[i] == 0)
      {
         continue;
      }
      else if (vt[i] == 2)
      {
         printf("method for nominal data will be implemented later\n");
         exit(1);
      }
      ll = gsl_matrix_get(JJ, i, 0) - 1 + 0.0000001;
      if (ll + 1 <= 2)
      {
         continue;
      }
      for (int j = 1; j < ll; j++)
      {
         a = gsl_matrix_get(Gama[i], j - 1, 0);
         b = gsl_matrix_get(Gama[i], j + 1, 0);

         for (int k = 0; k < N; k++)
         {
            yy = gsl_matrix_get(YY, k, i) + 0.0000001;
            z = gsl_matrix_get(ZZ, k, i);
            if ((yy == j) && (a < z))
            {
               a = z;
            }
            if ((yy == j + 1) && (b > z))
            {
               b = z;
            }
         }
         t1 = (b - a) * gsl_rng_uniform (r) + a;
         gsl_matrix_set(Gama[i], j, 0, t1);
      }
   }

   //free of memory
   gsl_matrix_free(ZZ);
}

// sample covariance and correlation matrix
// used by GS and GSM
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, const int & m0, const int & N)
{
   int K = Sigma -> size1;

   gsl_matrix ** WXB = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
   gsl_matrix ** T = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
   gsl_matrix * Txx = gsl_matrix_calloc(K, K);
   for (int i = 0; i < K; i++)
   {
      for (int j = 0; j < K; j++)
      {
         gsl_matrix_set(Txx, i, j, 0);
      }
   }
   for (int n = 0; n < N; n++)
   {
      WXB[n] = gsl_matrix_alloc(K, 1);
      T[n] = gsl_matrix_alloc(K, K);
   }
   double a, b, c;
   for (int n = 0; n < N; n++)
   {
      for (int k = 0; k < K; k++)
      {
         a = gsl_matrix_get(W, n * K + k, 0);
         b = gsl_matrix_get(XB, n * K + k, 0);
         c = a - b;
         gsl_matrix_set(WXB[n], k, 0, c);
      }

      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, WXB[n], WXB[n], 0.0, T[n]);

      gsl_matrix_add(Txx, T[n]);
   }
   gsl_matrix_add(Txx, Psigma);

   int dfinv = N + m0;
   InvWishart(Sigma, dfinv, Txx);
   CorrelationM(R, Sigma);

   gsl_matrix_free(Txx);
   for (int n = 0; n < N; n++)
   {
      gsl_matrix_free(WXB[n]);
      gsl_matrix_free(T[n]);
   }
   free(WXB);
   free(T);
}

// sample covariance and correlation matrix
// used by GS and GSM
void SampleSig_New(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, const int & m0, const int & N)
{
   int K = Sigma -> size1;

    gsl_matrix * WXB = gsl_matrix_alloc(K, 1);
    gsl_matrix * T = gsl_matrix_alloc(K, K);
    gsl_matrix * Txx = gsl_matrix_calloc(K, K);
    for (int i = 0; i < K; i++)
    {
      for (int j = 0; j < K; j++)
      {
         gsl_matrix_set(Txx, i, j, 0);
      }
   }
   double a, b, c;
   for (int n = 0; n < N; n++)
   {
      for (int k = 0; k < K; k++)
      {
         a = gsl_matrix_get(W, n * K + k, 0);
         b = gsl_matrix_get(XB, n * K + k, 0);
         c = a - b;
         gsl_matrix_set(WXB, k, 0, c);
      }

      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, WXB, WXB, 0.0, T);

      gsl_matrix_add(Txx, T);
   }

   gsl_matrix_add(Txx, Psigma);

   int dfinv = N + m0;
   InvWishart(Sigma, dfinv, Txx);
   CorrelationM(R, Sigma);

   gsl_matrix_free(Txx);
   gsl_matrix_free(WXB);
   gsl_matrix_free(T);
}
