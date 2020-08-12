// Parameter-extended random walk algorithm
// for longitudinal binary outcomes Y
// including missing values situation

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

gsl_rng * r;

void symmetric(gsl_matrix * A);
void Wishart(gsl_matrix * w, int m, gsl_matrix * sigma);
double rltruncnorm(double mu, double sd, double L);
double rrtruncnorm(double mu, double sd, double L);
double truncnorm(double mu, double sd, double a, double b);
void Inv(gsl_matrix * Ainverse, gsl_matrix * A);
void cholbeta(gsl_matrix * mu, gsl_matrix * cov, gsl_matrix * Z, gsl_matrix * X, gsl_matrix ** XX, gsl_matrix * Sigma, int N, gsl_matrix * b, gsl_matrix * C);
void multinorm(gsl_matrix * Mnorm, gsl_matrix * mu, gsl_matrix * sigma);
void Mat_I_J(gsl_matrix * B, gsl_matrix * A, int k);
void MatI_J(gsl_matrix * B, gsl_matrix * A, int k);
void Mat_IJ(gsl_matrix * B, gsl_matrix * A, int k);
void Mat_I(gsl_matrix * B, gsl_matrix * A, int k);
void SampleZ(gsl_matrix * Z, gsl_matrix * sigma, gsl_matrix * XB, gsl_matrix * Y, int N, gsl_matrix ** Gama, gsl_matrix * JJ);
void SampleGama(int N, int rep, gsl_matrix * YY, gsl_matrix * Z, gsl_matrix **Gama, gsl_matrix * JJ);
void CorrelationM(gsl_matrix * R, gsl_matrix * W);
double Deter(gsl_matrix * A);
double Priornew(int m0, gsl_matrix * sigma, gsl_matrix * W);
double Posterior(int N, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB);
void InvWishart(gsl_matrix * SigmaMat, int m, gsl_matrix * ScaleMat);
int MHsampler(gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma, int N, int m, int m0);
void Gibbsampler(gsl_matrix * Y, gsl_matrix * X, gsl_matrix * YY, gsl_matrix ** XX, gsl_matrix * Psigma, gsl_matrix ** Gama, gsl_matrix * JJ, int N, int K, int rep, int Inum, int m, int m0,gsl_matrix * b, gsl_matrix * C,
				 char * outBetaFile, char * outZFile, char * outRFile, char * outWRFile, char * outGamaFile);

int main(void)
{
	int i, j, k;
	double a;
	//regarding the prior of correlation R
	//double rh0 = 0.5;

	// N: total number of individuals
	// K: number of regression parameters
	// m0: degrees of freedom of prior for correlation
	// m: degrees of freedom of proposed density of correlation
	// Inum: number of iterations
	// rep: number of repeated measures for each individual
    // This part will be changed later, all parameters will be input from program
	int N = 50;
	int K = 2;
	int rep = 5;
	int m0 = 10;
	int m = 100;
	int Inum = 10000;
	int ordcut = 4;

	// X: covariate matrix which is N*r by K
	// Y: observation vector which is treated as a N*r by 1 matrix
	// **Z: latent vectors which is also treated as a N*r by 1 matrix
	// XB: product of X and regression parameters beta
	// **beta: gsl_matrix which is regression parameters
	// **R: gsl_matrix which is correlation matrix of observed Y for each individual
	// b and C are the prior of beta
	gsl_matrix * XB, * X, * Y, * YY, * Psigma, ** XX, ** Gama, * b , * C;

	// allocate memory
	XB = gsl_matrix_alloc(N * rep, 1);
    X = gsl_matrix_alloc(N * rep, K);
	Y = gsl_matrix_alloc(N * rep, 1);
	YY = gsl_matrix_alloc(N, rep);
	Psigma = gsl_matrix_alloc(rep, rep);
	b = gsl_matrix_alloc(K, 1);
	C = gsl_matrix_calloc(K, K);
	XX = (gsl_matrix **) calloc(N, sizeof(gsl_matrix * ));

	for(k = 0; k < N; k++)
	{
		XX[k] = gsl_matrix_alloc(rep, K);
	}

	gsl_matrix *JJ = gsl_matrix_alloc(rep,1);
    for(k=0; k<rep; k++)
    {
        gsl_matrix_set(JJ, k, 0, ordcut);
    }
   	//Gama[k] means all the cutpoints for kth measures
    Gama = (gsl_matrix **) calloc(rep, sizeof(gsl_matrix*));
	for(k = 0; k < rep; k++)
	{
		Gama[k] = gsl_matrix_alloc(gsl_matrix_get(JJ, k, 0),1);
	}
	//set the sigma matrix for the prior distribution of R and D.
    //This is to set the prior of R to be the identity matrix
	gsl_matrix_set_identity(Psigma);
	gsl_matrix_scale(Psigma, m0);
	//Set the prior for beta

	for(i = 0; i < K; i++)
	{
        gsl_matrix_set(b, i, 0, 0);
    }
    gsl_matrix_set_zero (C);
    //gsl_matrix_set_identity(C);
    //gsl_matrix_scale(C, 10000);
    // set the seed for random number generator
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

    // input data X, Y from files
    char  *dirName, * fileName;
    FILE * fp;
    dirName = new char [80];
    fileName = new char [80];
    strcpy(dirName, "../DataResults/");
    strcpy(fileName, dirName);
    strcat(fileName, "outY.dat");
 	if((fp = fopen(fileName, "r")) == NULL)
	{
		printf("can not open file %s!\n", fileName);
		exit(1);
	}
	gsl_matrix_fscanf(fp, Y);
    fclose(fp);
    for(i=0;i<N;i++)
	{
		for(j=0;j<rep;j++)
		{
			a = gsl_matrix_get(Y, i*rep+j, 0);
			gsl_matrix_set(YY, i, j, a);
		}
	}
	strcpy(fileName, dirName);
    strcat(fileName, "outX.dat");
 	if((fp = fopen(fileName, "r")) == NULL)
	{
		printf("can not open file %s!\n", fileName);
		exit(1);
	}
	gsl_matrix_fscanf(fp, X);
    fclose(fp);
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < rep; j++)
		{
			for(k = 0; k < K; k++)
			{
				a = gsl_matrix_get(X, i * rep + j, k);
				gsl_matrix_set(XX[i], j, k, a);
			}
		}
	}
    // set up output files`
	char * outBetaFile, * outZFile, * outRFile, * outWRFile, * outGamaFile;
	outBetaFile = new char [80];
    outZFile = new char [80];
    outRFile = new char [80];
    outWRFile = new char [80];
    outGamaFile = new char [80];
	strcpy(outBetaFile, dirName);
	strcpy(outZFile, dirName);
	strcpy(outRFile, dirName);
	strcpy(outWRFile, dirName);
	strcpy(outGamaFile, dirName);
	strcat(outBetaFile, "CorrOrdBeta.dat");
	strcat(outZFile, "CorrOrdZ.dat");
	strcat(outRFile, "CorrOrdR.dat");
	strcat(outWRFile, "CorrOrdWR.dat");
	strcat(outGamaFile, "CorrOrdGama.dat");
	// Do Gibbs Sampling based on the above data.
	Gibbsampler(Y, X, YY, XX, Psigma, Gama, JJ, N, K, rep, Inum, m, m0, b, C, outBetaFile, outZFile, outRFile, outWRFile, outGamaFile);

	// free memory;
	gsl_matrix_free(X);
	gsl_matrix_free(Y);
	gsl_matrix_free(YY);
	gsl_matrix_free(XB);
	gsl_matrix_free(Psigma);
	gsl_matrix_free(b);
	gsl_matrix_free(C);
	gsl_matrix_free(JJ);
	for(k =0; k < N; k++)
	{
		gsl_matrix_free(XX[k]);
	}
	free(XX);
	for(k=0;k<rep;k++)
	{
		gsl_matrix_free(Gama[k]);
	}
	free(Gama);
	delete []fileName;
	delete []dirName;
	delete []outBetaFile;
	delete []outZFile;
	delete []outRFile;
	delete []outWRFile;
	delete []outGamaFile;
}
// make a matrix to be symmetric
// due to numerical computation, a symmetric matrix may become non-symmetric
void symmetric(gsl_matrix * A)
{
	int i, j, k1, k2;
	double a, b, d;

	k1= A -> size1;
	k2= A -> size2;

	if(k1 != k2)
	{
		printf("the matrix is not square matrix!\n");
		exit(1);
	}

	for(i = 0; i< k1; i++)
	{
		for(j = 0; j < i; j++)
		{
			a = gsl_matrix_get(A, j, i);
			b = gsl_matrix_get(A, i, j);
			d = (a + b ) / 2;
			gsl_matrix_set(A, i, j, d);
			gsl_matrix_set(A, j, i, d);
		}
	}
}
// get the inverse of the matrix A which is denoted as Ainverse without changing A
void Inv(gsl_matrix * Ainverse, gsl_matrix * A)
{
	int signum;
    int k1, k2;
	k1 = A -> size1;
	k2 = A -> size2;
	if(k1 != k2)
	{
		printf("the matrix is not a square matrix!\n");
		exit(1);
	}
	gsl_permutation * p = gsl_permutation_alloc(k1);
	gsl_matrix * Q = gsl_matrix_alloc(k1, k1);
	gsl_matrix_memcpy(Q, A);
    gsl_linalg_LU_decomp(Q, p, &signum);
	gsl_linalg_LU_invert(Q, p, Ainverse);

	// free of memory;
	gsl_matrix_free(Q);
	gsl_permutation_free(p);
}
// generate a random matrix from Wishart Distribution with df = m and sigma.
// all the matrix have been allocated and sigma will not be changed
void Wishart(gsl_matrix * w, int m, gsl_matrix * sigma)
{
    int k;
    k = sigma -> size1;
    gsl_matrix * Q = gsl_matrix_alloc(k, k);
    gsl_matrix *work = gsl_matrix_alloc(k, k);

    gsl_matrix_memcpy(Q, sigma);
    gsl_linalg_cholesky_decomp(Q);
    gsl_ran_wishart(r, m, Q, w, work);

    gsl_matrix_free(Q);
    gsl_matrix_free(work);
}
//get correlation matrix R from W
void CorrelationM(gsl_matrix * R, gsl_matrix * W)
{
	int i, j, k;
	double a, b, c;
	k = W -> size1;

	for(i = 0; i < k; i++)
	{
		for(j = 0; j < k; j++)
		{
			a = gsl_matrix_get(W, i, i);
			b = gsl_matrix_get(W, j, j);
			c = gsl_matrix_get(W, i, j);
			gsl_matrix_set(R, i, j, c / sqrt(a * b));
		}
	}
}
// generate random variable from left truncated normal distribution with mean u and standard
// variance std at value of L
double rltruncnorm(double mu, double sd, double L)
{
	double t, x, y;
	t = 1 - gsl_sf_erf_Q((L - mu) / sd);
    x = (1 - t) * gsl_rng_uniform(r) + t;
	y = gsl_cdf_ugaussian_Pinv(x) * sd + mu;
	return(y);
}
//generate random variable from right truncated normal distribution with mean u and
// standard variance std at value of L
double rrtruncnorm(double mu, double sd, double L)
{
	double t, x, y;
	t = 1 - gsl_sf_erf_Q((L - mu) / sd);
	x = t * gsl_rng_uniform(r);
	y = gsl_cdf_ugaussian_Pinv(x) * sd + mu;
	return(y);
}
//generate random variable from left and right truncated normal distribution with mean u
//standard variance std at value of [a,b]
double truncnorm(double mu, double sd, double a, double b)
{
	double t1, t2, x, y;
	t1= 1- gsl_sf_erf_Q ((a-mu)/sd);
	t2= 1- gsl_sf_erf_Q ((b-mu)/sd);
	x = (t2-t1)*gsl_rng_uniform(r) + t1;
	y = gsl_cdf_ugaussian_Pinv(x) * sd + mu;
	return(y);
}
// get the mean vector of "mu" and co-variance matrix "cov" of beta
// Input latent variable Z, covariate X, sigma, and sample size N
void cholbeta(gsl_matrix * mu, gsl_matrix * cov, gsl_matrix * Z, gsl_matrix * X, gsl_matrix ** XX, gsl_matrix * Sigma, int N, gsl_matrix * b, gsl_matrix * C)
{
	int k, n2, i, j;
	double a;
	symmetric(Sigma);
	k = Sigma -> size1;
	n2 = X -> size2;
	gsl_matrix ** ZZ;
	ZZ = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
	for(i = 0; i < N; i++)
	{
        ZZ[i] = gsl_matrix_alloc(k, 1);
	}
	gsl_matrix * Q = gsl_matrix_alloc(k, k);
	gsl_matrix * Sig = gsl_matrix_calloc(k, k);
	gsl_matrix * Siginv = gsl_matrix_alloc(k, k);
	gsl_matrix * T1 = gsl_matrix_alloc(n2, k);
	gsl_matrix * T2 = gsl_matrix_alloc(n2, n2);
	gsl_matrix * T3 = gsl_matrix_alloc(k, 1);
	gsl_matrix * T4 = gsl_matrix_alloc(n2, 1);
	gsl_matrix * Txx = gsl_matrix_calloc(n2, n2);
	gsl_matrix * Txz = gsl_matrix_calloc(n2, 1);
	gsl_matrix * InvC = gsl_matrix_calloc(n2, n2);
	gsl_matrix * CIb = gsl_matrix_calloc(n2, 1);

	// copy sigma to Q, so sigma will not be changed;
	gsl_matrix_memcpy(Q, Sigma);

	// choleskey decomposition
	gsl_linalg_cholesky_decomp(Q);

	// get the lower triangle matrix Sig
	for(i = 0; i < k; i++)
	{
	  for(j = 0; j <= i; j++)
	  {
		  a = gsl_matrix_get(Q, i, j);
		  gsl_matrix_set(Sig, i, j, a);
	  }
	}
	Inv(Siginv, Sig);

	for(i = 0;i < N; i++)
	{
		for(j = 0; j < k; j++)
		{
			a = gsl_matrix_get(Z, i * k + j,0);
			gsl_matrix_set(ZZ[i], j, 0, a);
		}

		gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, XX[i], Siginv, 0.0, T1);

		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T1, T1, 0.0, T2);

		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Siginv, ZZ[i], 0.0, T3);

		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T1, T3, 0.0, T4);

		gsl_matrix_add(Txx, T2);

		gsl_matrix_add(Txz, T4);
	}

	a = gsl_matrix_get(C, 0, 0);

    if(a==0) gsl_matrix_memcpy(InvC, C);

    else if (a>0) Inv(InvC, C);

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
	gsl_matrix_free(InvC);
	gsl_matrix_free(CIb);
	for(i = 0; i < N; i++)
	{
		gsl_matrix_free(ZZ[i]);
	}
	free(ZZ);
}
// get the random number from multivariate normal with mean mu and covariance matrix sigma
// mu and sigma will not be changed in this routine
// In this routine, the random vector is treated as a matrix.
void multinorm(gsl_matrix * Mnorm, gsl_matrix * mu, gsl_matrix * sigma)
{
	int i, k;
	k = sigma -> size1;

	//  Cholesky factor of sigma
    gsl_matrix * Q = gsl_matrix_alloc(k, k);
	symmetric(sigma);
	gsl_matrix_memcpy(Q, sigma);
    gsl_linalg_cholesky_decomp(Q);

    gsl_vector * mymu = gsl_vector_alloc(k);
    gsl_vector * res = gsl_vector_alloc(k);
    double a;
    for(i = 0; i < k; i++)
    {
        a = gsl_matrix_get(mu, i, 0);
        gsl_vector_set(mymu, i, a);
    }

    gsl_ran_multivariate_gaussian(r, mymu, Q, res);

    for(i = 0; i < k; i++)
    {
        a = gsl_vector_get(res, i);
        gsl_matrix_set(Mnorm, i, 0, a);
    }
	// free memory
	gsl_matrix_free(Q);
	gsl_vector_free(mymu);
	gsl_vector_free(res);
}
// get the matrix B which is without kth row and kth column of A.
void Mat_I_J(gsl_matrix * B, gsl_matrix * A, int k)
{
	int i, j;
	int m, n;
	m = A -> size1;
	n = A -> size2;
	for(i = 0; i < m - 1; i++)
	{
		for(j = 0; j < n - 1; j++)
		{
			if((i < k) && (j < k))
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i, j));
			else if((i < k) && (j >= k))
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i, j + 1));
			else if((i >= k) && (j < k))
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i + 1, j));
			else if((i >= k) && (j >= k))
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i + 1, j + 1));
		}
	}
}
// get the kth row from matrix A with kth element deleted.
void MatI_J(gsl_matrix * B, gsl_matrix * A, int k)
{
	int i, m;
	m = A -> size2;
	for(i = 0; i < m - 1; i++)
	{
		if(i < k)
			gsl_matrix_set(B, 0, i, gsl_matrix_get(A, k, i));
		else
			gsl_matrix_set(B, 0, i, gsl_matrix_get(A, k, i + 1));
	}
}
//get the kth column from matrix A with kth element deleted.
void Mat_IJ(gsl_matrix * B, gsl_matrix * A, int k)
{
	int i, m;
	m = A -> size1;
	for(i = 0; i < m - 1; i++)
	{
		if(i < k)
			gsl_matrix_set(B, i, 0, gsl_matrix_get(A, i, k));
		else
			gsl_matrix_set(B, i, 0, gsl_matrix_get(A, i + 1, k));
	}
}
// get the matrix B from matrix A with the kth row deleted.
void Mat_I(gsl_matrix * B, gsl_matrix * A, int k)
{
	int i, j, m, n;
	m = A -> size1;
	n = A -> size2;
	for(i = 0; i < m - 1; i++)
	{
		for(j = 0; j < n; j++)
		{
			if(i < k)
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i, j));
			else
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i + 1, j));
		}
	}
}
// get the sample of latent variable Z
// Input:
//	sigma: m*m covariance matrix;
//	XB: (N*m)*1 matrix;
//  Y: observed discrete outcomes;
//  N: sample size.
void SampleZ(gsl_matrix * Z, gsl_matrix * sigma, gsl_matrix * XB, gsl_matrix * Y, int N, gsl_matrix ** Gama, gsl_matrix * JJ)
{
	int i, j, k, n, m;
	double sigma22, sigma221, u00, a, b, num;
	m = sigma -> size1;
	gsl_matrix * sigma11 = gsl_matrix_alloc (m-1, m-1);
	gsl_matrix * Isigma11 = gsl_matrix_alloc (m-1, m-1);
	gsl_matrix * sigma21 = gsl_matrix_alloc (m-1,1);
	gsl_matrix * T1 = gsl_matrix_alloc (1,m-1);
	gsl_matrix * T2 = gsl_matrix_alloc (1,1);
	gsl_matrix * T3 = gsl_matrix_alloc (m-1,1);
	gsl_matrix * T4 = gsl_matrix_alloc (m-1,1);
	gsl_matrix * T5 = gsl_matrix_alloc (1,1);

	for(i=0; i<N*m; i++)
	{
		//aa = gsl_matrix_get(Y,i,0);
		j = i % m ;
		k = (i-j)/m;
		num = gsl_matrix_get(JJ, j, 0);
		sigma22 = gsl_matrix_get(sigma,j,j);
		Mat_I_J(sigma11, sigma, j);
		Inv(Isigma11, sigma11);
		Mat_IJ(sigma21, sigma, j);

		gsl_blas_dgemm (CblasTrans,CblasNoTrans,1.0, sigma21, Isigma11,0.0, T1);
		gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1.0, T1, sigma21,0.0, T2);
		b = gsl_matrix_get(T2,0,0);
		sigma221 = sqrt(sigma22 - b);
		for(n=0; n<m-1; n++)
		{
				if((k*m+n)<i)
				{	a = gsl_matrix_get(Z, k*m+n, 0);
					gsl_matrix_set(T3, n, 0, a);
				}
				else
				{
					a = gsl_matrix_get(Z, k*m+n+1, 0);
					gsl_matrix_set(T3, n, 0, a);
				}
		}
		for(n=0; n<m-1; n++)
		{
				if((k*m+n)<i)
				{
					a = gsl_matrix_get(XB, k*m+n, 0);
					gsl_matrix_set(T4, n, 0, a);
				}
				else
				{
					a = gsl_matrix_get(XB, k*m+n+1, 0);
					gsl_matrix_set(T4, n, 0, a);
				}
		}

		gsl_matrix_sub (T3, T4);

		gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1.0, T1,T3,0.0, T5);

		u00 = gsl_matrix_get(XB, i, 0) + gsl_matrix_get(T5, 0, 0);


		if(gsl_matrix_get(Y ,i, 0) == 0)
		{
				a = rrtruncnorm(u00, sigma221, 0.0);
				gsl_matrix_set(Z, i, 0, a);
		}

		else if (gsl_matrix_get(Y, i, 0)==num-1)
		{
				a = rltruncnorm(u00, sigma221, gsl_matrix_get(Gama[j], num-2, 0));
				gsl_matrix_set(Z, i, 0, a);
		}
		else if (gsl_matrix_get(Y, i, 0)==999)
		{
				a = gsl_ran_gaussian (r, sigma221) + u00;
				gsl_matrix_set(Z, i, 0, a);
		}
		else
		{
				for(n=1; n<num-1; n++)
				{
					if(gsl_matrix_get(Y ,i, 0)==n)
					{
						a = truncnorm(u00, sigma221, gsl_matrix_get(Gama[j], n-1, 0), gsl_matrix_get(Gama[j], n,0));
						gsl_matrix_set(Z, i, 0, a);
					}

				}
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
//sample the Gama matrix
void SampleGama(int N, int rep, gsl_matrix * YY, gsl_matrix * Z, gsl_matrix **Gama, gsl_matrix * JJ)
{
	int i,j,k;
	double a, b, t1, y, z;
    gsl_matrix * ZZ = gsl_matrix_alloc(N,rep);
	for(i=0;i<N;i++)
	{
		for(j=0;j<rep;j++)
		{
			a = gsl_matrix_get(Z, i*rep+j, 0);
			gsl_matrix_set(ZZ, i, j, a);
		}
	}
	for(i=0; i<rep; i++)
	{
		for(j=1; j<gsl_matrix_get(JJ,i,0)-1; j++)
		{
			a = gsl_matrix_get(Gama[i], j-1,0);
			b = gsl_matrix_get(Gama[i], j+1,0);

			for(k=0;k<N;k++)
			{
				y = gsl_matrix_get(YY, k, i);
				z = gsl_matrix_get(ZZ, k, i);
				if((y==j)&&(a<z))
				{
					a=z;
				}
				if((y==j+1)&&(b>z))
				{
					b=z;
				}
			}
			t1 = (b-a) * gsl_rng_uniform (r) + a;
			gsl_matrix_set(Gama[i], j, 0, t1);
		}
	}
	//free of memory
	gsl_matrix_free(ZZ);
}
//computing the log determinant of a positive definite matrix A.
//In this routine, A is not changed.
double Deter(gsl_matrix * A)
{
	int k, signum;
	k = A -> size1;
	gsl_matrix * Q = gsl_matrix_alloc(k, k);
    gsl_permutation * p = gsl_permutation_alloc(k);

	symmetric(A);
	gsl_matrix_memcpy(Q, A);
    gsl_linalg_LU_decomp(Q, p, &signum);
    double sum;
    sum = log(gsl_linalg_LU_det(Q, signum));

	gsl_matrix_free(Q);
	gsl_permutation_free(p);

	return(sum);
}
double Priornew(int m0, gsl_matrix * sigma, gsl_matrix * W)
{
	int k, i;
	double a1,a2,b1, b2, t1, t2;
	k = W -> size1;

	gsl_matrix * InvW = gsl_matrix_alloc(k, k);

    gsl_matrix * A = gsl_matrix_alloc(k, k);

	a1 = Deter(W);
	a2 = Deter(sigma);

	Inv(InvW, W);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sigma, InvW, 0.0, A);

	t1 = 0;
    t2 = 0;
	for(i = 0; i < k; i++)
	{
		b1 = gsl_matrix_get(A, i, i);
		b2 = gsl_matrix_get(W, i, i);
		t1 = t1 + b1;
		t2 = t2 + log(b2);
    }

	//free memory
	//gsl_matrix_free(Invsigma);
	gsl_matrix_free(InvW);
	gsl_matrix_free(A);
	//compute the Jacobian
	return(0.5 * (k - 1) * t2 - 0.5 * (m0 + k + 1) * a1 - 0.5 * t1 + 0.5 * m0 * a2);
}
//computing the log posterior density of the correlation matrix not including the prior part.
double Posterior(int N, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB)
{
	int i, j, k;
	double a, sum1, sum2, Txx = 0.0;
	k = R -> size1;
	gsl_matrix * Q1 = gsl_matrix_alloc(k * N, 1);
	gsl_matrix * Q2 = gsl_matrix_alloc(k, k);
	gsl_matrix * Sig = gsl_matrix_calloc(k, k);
	gsl_matrix * Siginv = gsl_matrix_alloc(k, k);
	gsl_matrix ** QQ;
	QQ = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
	for(i = 0; i < N; i++)
	{
		QQ[i] = gsl_matrix_alloc(k, 1);
	}

	gsl_matrix * T1 = gsl_matrix_alloc(1, k);
	gsl_matrix * T2 = gsl_matrix_alloc(1, 1);

	gsl_matrix_memcpy(Q1, Z);
	gsl_matrix_sub(Q1, XB);

	gsl_matrix_memcpy(Q2, R);

	// choleskey decomposition
	gsl_linalg_cholesky_decomp(Q2);

	// get the lower triangle matrix Sig
	for(i = 0; i < k; i++)
	{
	  for(j = 0; j <= i; j++)
	  {
		  gsl_matrix_set(Sig, i, j, gsl_matrix_get(Q2, i, j));
	  }
	}

	Inv(Siginv, Sig);
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < k; j++)
		{
			a = gsl_matrix_get(Q1, i * k + j,0);
			gsl_matrix_set(QQ[i], j, 0, a);
		}

		gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, QQ[i], Siginv, 0.0, T1);

		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T1, T1, 0.0, T2);

		Txx = Txx + gsl_matrix_get(T2, 0, 0);
	}

	sum1 = (-N/2) * Deter(R);
	sum2 = (-1.0 / 2.0) * Txx;

	// free memory;
	gsl_matrix_free(Q1);
	gsl_matrix_free(Q2);
	gsl_matrix_free(Sig);
	gsl_matrix_free(Siginv);
	gsl_matrix_free(T1);
	gsl_matrix_free(T2);
	for(i = 0; i < N; i++)
	{
		gsl_matrix_free(QQ[i]);
	}
	free(QQ);
	return(sum1 + sum2);
}
// all the matrix have been allocated and sigma will not be changed
void InvWishart(gsl_matrix * SigmaMat, int m, gsl_matrix * ScaleMat)
{
    int k;
    k = ScaleMat -> size1;
    gsl_matrix * Q = gsl_matrix_alloc(k, k);
    gsl_matrix *work = gsl_matrix_alloc(k, k);
    gsl_matrix *InvScale = gsl_matrix_alloc(k, k);
    gsl_matrix *InvSig = gsl_matrix_alloc(k, k);

    Inv(InvScale, ScaleMat);
    int df;
    //df = m-k-1;
    df = m;
    gsl_matrix_memcpy(Q, InvScale);
    gsl_linalg_cholesky_decomp(Q);
    gsl_ran_wishart(r, df, Q, InvSig, work);

    Inv(SigmaMat, InvSig);

    gsl_matrix_free(Q);
    gsl_matrix_free(work);
    gsl_matrix_free(InvScale);
    gsl_matrix_free(InvSig);
}
//M-H sampling algorithm.
//The information for input is described in the main function.
int MHsampler(gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma, int N, int m, int m0)
{
	int k, df;
	k = R -> size1;
	double u, p1, p2, post1, post2, pro1, pro2;
	double rho;
	gsl_matrix * Wc = gsl_matrix_alloc(k, k);
	gsl_matrix * Rc = gsl_matrix_alloc(k, k);
	gsl_matrix * Wcc = gsl_matrix_alloc(k, k);
	gsl_matrix * Wrr = gsl_matrix_alloc(k, k);

	df=m-k-1;
	gsl_matrix_memcpy(Wrr, WR);
	gsl_matrix_scale(Wrr, df);

    //Wishart(Wc, m, Wrr);
    InvWishart(Wc, m, Wrr);
	gsl_matrix_memcpy(Wcc, Wc);
	gsl_matrix_scale(Wcc, df);

	CorrelationM(Rc, Wc);

    u = gsl_rng_uniform(r);

	p1 = Priornew(m0, Psigma, Wc);
	p2 = Priornew(m0, Psigma, WR);
	post1 = Posterior(N, Rc, Z, XB);
	post2 = Posterior(N, R, Z, XB);
	pro1 = Priornew(m, Wcc, WR);
	pro2 = Priornew(m, Wrr, Wc);

    rho = p1 + post1 + pro1 - p2 - post2 - pro2;

	if(u <= exp(rho))
	{
		gsl_matrix_memcpy(WR, Wc);
		gsl_matrix_memcpy(R, Rc);
		gsl_matrix_free(Wc);
		gsl_matrix_free(Rc);
		gsl_matrix_free(Wcc);
		gsl_matrix_free(Wrr);
		return(1);
	}
	gsl_matrix_free(Wc);
	gsl_matrix_free(Rc);
	gsl_matrix_free(Wcc);
	gsl_matrix_free(Wrr);
	return(0);
}
// Gibbs Sampling algorithm
// Input:
//	Y: (N * K) * 1 matrix
//	X: (N * K)* n2 matrix
//	N: the number of individuals
//	K: the number of parameters
//  r: the number of repeated measures
//	Inum: number of total of iterations
//	m0: prior parameter for correlation matrix R
//	m: proposed parameter for sampled correlation matrix R
//	lamda: jumping parameter
//	epsilon: precision
//	beta: an array of regression parameter for each iteration, each element is n2 * 1 matrix
//	R: an array of correlation matrix for each iteration, each element is K * K matrix
//	Z: an array of latent variable for each iteration, each element is (N * K) * 1 matrix
void Gibbsampler(gsl_matrix * Y, gsl_matrix * X, gsl_matrix * YY, gsl_matrix ** XX, gsl_matrix * Psigma, gsl_matrix ** Gama, gsl_matrix * JJ, int N, int K, int rep, int Inum, int m, int m0,gsl_matrix * b, gsl_matrix * C,
				 char * outBetaFile, char * outZFile, char * outRFile, char * outWRFile, char * outGamaFile)
{
	int i, k, sum = 0;
	double a;
	gsl_matrix * mu = gsl_matrix_alloc(K, 1);
	gsl_matrix * cov = gsl_matrix_alloc(K, K);
	gsl_matrix * XB = gsl_matrix_alloc(N * rep, 1);

	// here we define Z, beta and R to replace the previous parameters;
	gsl_matrix *newZ = gsl_matrix_alloc(N * rep, 1);
	gsl_matrix *newR = gsl_matrix_alloc(rep, rep);
	gsl_matrix *newBeta = gsl_matrix_alloc(K, 1);
	gsl_matrix *newWR = gsl_matrix_alloc(rep, rep);
    gsl_matrix *RR = gsl_matrix_alloc(rep, rep);

	// set the initial value for R;
	//gsl_matrix_set_identity(RR);
	gsl_matrix_set_identity(newR);
	gsl_matrix_set_identity(newWR);

	gsl_matrix **newGama;
	//initialize newGama
	newGama = (gsl_matrix **) calloc(rep, sizeof(gsl_matrix*));
	for(k=0; k<rep; k++)
	{
		newGama[k] = gsl_matrix_calloc(gsl_matrix_get(JJ,k,0),1);
		if(gsl_matrix_get(JJ,k,0)>1)
			gsl_matrix_set(newGama[k], gsl_matrix_get(JJ,k,0)-1, 0, 10000.0);
	}

	// set the initial value for Z;
    for(i = 0; i < N * rep; i++)
	{
		a = gsl_matrix_get(Y, i, 0);
		if(a > 0)
			gsl_matrix_set(newZ, i, 0, rltruncnorm(0.0, 1.0, 0.0));
		//else if(fabs(a-999)<0.01)
		//	gsl_matrix_set(newZ, i, 0, gsl_ran_gaussian(r, 1.0));
		else
			gsl_matrix_set(newZ, i, 0, rrtruncnorm(0.0, 1.0, 0.0));
	}

	// open the file for output;
	FILE *betafp, *zfp, *rfp, *wrfp, *gamafp;
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
	if((rfp = fopen(outRFile,"w")) == NULL)
	{
		printf("can not open file %s in Gibbssampler()!\n", outRFile);
		exit(1);
	}
	if((wrfp = fopen(outWRFile,"w")) == NULL)
	{
		printf("can not open file %s in Gibbssampler()!\n", outWRFile);
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

	for(i = 1; i <= Inum; i++)
	{
	   if((i % 100) == 0)
       {
	     printf("step in Gibbs Sampling : %d\n", i);
       }
   	   cholbeta(mu, cov, newZ, X, XX, newR, N, b, C);
   	   multinorm(newBeta, mu, cov);
       gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, newBeta, 0.0, XB);
       sum = sum + MHsampler(newWR, newR, newZ, XB, Psigma, N, m, m0);
       SampleZ(newZ, newR, XB, Y, N, newGama, JJ);
       SampleGama(N, rep, YY, newZ, newGama, JJ);

		// this part is for the output;
        gsl_matrix_fprintf(betafp, newBeta, "%g");
        //if((i % 100) == 0) gsl_matrix_fprintf(zfp, newZ, "%g");
        gsl_matrix_fprintf(rfp, newR, "%g");
        gsl_matrix_fprintf(wrfp, newWR, "%g");
        //gsl_matrix_fprintf(gamafp, newGama, "%g");
        for(k=0; k<rep; k++)
        {
               gsl_matrix_fprintf(gamafp, newGama[k], "%g");
        }
    }
	printf("%d\n", sum);
	// close the file;
	fclose(betafp);
	fclose(zfp);
	fclose(rfp);
	fclose(wrfp);
	fclose(gamafp);
	// free memory;
	gsl_matrix_free(RR);
	gsl_matrix_free(mu);
	gsl_matrix_free(cov);
	gsl_matrix_free(XB);
	gsl_matrix_free(newZ);
	gsl_matrix_free(newR);
	gsl_matrix_free(newBeta);
	gsl_matrix_free(newWR);
	for(k=0;k<rep;k++)
	{
		gsl_matrix_free(newGama[k]);
	}
	free(newGama);
}
