/////////////////////////////////////
// Name: Random.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related with random variables

#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include "Matrix.h"

extern gsl_rng * r;

// generate random variable from left truncated normal distribution with mean mu 
// and standard deviation sd at value of L
double rltruncnorm(const double & mu, const double & sd, const double & L)
{
	double t = 1 - gsl_sf_erf_Q((L - mu) / sd);
    double x = (1 - t) * gsl_rng_uniform(r) + t;
	double y = gsl_cdf_ugaussian_Pinv(x) * sd + mu;
	return(y);
}

// generate random variable from right truncated normal distribution with mean mu
// and standard deviation sd at value of L
double rrtruncnorm(const double & mu, const double & sd, const double & L)
{
	double t = 1 - gsl_sf_erf_Q((L - mu) / sd);
	double x = t * gsl_rng_uniform(r);
	double y = gsl_cdf_ugaussian_Pinv(x) * sd + mu;
	return(y);
}

// generate random variable from left and right truncated normal distribution with mean mu
// standard deviation sd at value of [a,b]
double truncnorm(const double & mu, const double & sd, const double & a, const double & b)
{
	double t1 = 1 - gsl_sf_erf_Q ((a - mu)/sd);
	double t2 = 1 - gsl_sf_erf_Q ((b - mu)/sd);
	double x = (t2 - t1) * gsl_rng_uniform(r) + t1;
	double y = gsl_cdf_ugaussian_Pinv(x) * sd + mu;
	return(y);
}

// generate a random matrix from inv Wishart Distribution with d.f. = m and sigma.
// all the matrix have been allocated and sigma will not be changed
void InvWishart(gsl_matrix * SigmaMat, const int & m, gsl_matrix * ScaleMat)
{
    int k = ScaleMat -> size1;
    gsl_matrix * Q = gsl_matrix_alloc(k, k);
    gsl_matrix * work = gsl_matrix_alloc(k, k);
    gsl_matrix * InvScale = gsl_matrix_alloc(k, k);
    gsl_matrix * InvSig = gsl_matrix_alloc(k, k);

    Inv(InvScale, ScaleMat);
    int df = m - k - 1;
    gsl_matrix_memcpy(Q, InvScale);
    gsl_linalg_cholesky_decomp(Q);
    gsl_ran_wishart(r, df, Q, InvSig, work);
    Inv(SigmaMat, InvSig);

    gsl_matrix_free(Q);
    gsl_matrix_free(work);
    gsl_matrix_free(InvScale);
    gsl_matrix_free(InvSig);
}

// generate a random matrix from Wishart Distribution with d.f. = m and sigma.
// all the matrix have been allocated and sigma will not be changed
void Wishart(gsl_matrix * w, const int & m, gsl_matrix * Sigma)
{
    int k = Sigma -> size1;
    gsl_matrix * Q = gsl_matrix_alloc(k, k);
    gsl_matrix * work = gsl_matrix_alloc(k, k);

    gsl_matrix_memcpy(Q, Sigma);
    gsl_linalg_cholesky_decomp(Q);
    gsl_ran_wishart(r, m, Q, w, work);

    gsl_matrix_free(Q);
    gsl_matrix_free(work);
}

// generate random number from multivariate normal with mean mu and covariance matrix sigma
// mu and sigma will not be changed in this function
// In this function, the random vector is treated as a matrix.
void multinorm(gsl_matrix * Mnorm, gsl_matrix * mu, gsl_matrix * sigma)
{
	int k = sigma -> size1;

	//  Cholesky factor of sigma
    gsl_matrix * Q = gsl_matrix_alloc(k, k);
	symmetric(sigma);
	gsl_matrix_memcpy(Q, sigma);
    gsl_linalg_cholesky_decomp(Q);

    gsl_vector * mymu = gsl_vector_alloc(k);
    gsl_vector * res = gsl_vector_alloc(k);
    double a;
    for(int i = 0; i < k; i++)
    {
        a = gsl_matrix_get(mu, i, 0);
        gsl_vector_set(mymu, i, a);
    }

    gsl_ran_multivariate_gaussian(r, mymu, Q, res);

    for(int i = 0; i < k; i++)
    {
        a = gsl_vector_get(res, i);
        gsl_matrix_set(Mnorm, i, 0, a);
    }

	// free memory
	gsl_matrix_free(Q);
	gsl_vector_free(mymu);
	gsl_vector_free(res);
}
