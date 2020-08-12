// Parameter Expanded Algorithm
// Standard Gibbs Sampling
// Multivariate Binary Outcomes Y
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
void InvWishart(gsl_matrix * SigmaMat, int m, gsl_matrix * ScaleMat);
void CorrelationM(gsl_matrix * R, gsl_matrix * Sigma);
void SampleSD(gsl_matrix * SD, gsl_matrix * R, gsl_matrix * Psigma, int m);
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, int m, int N);
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
void Gibbsampler(gsl_matrix * Y, gsl_matrix * X, gsl_matrix * YY, gsl_matrix ** XX, gsl_matrix * Psigma, gsl_matrix ** Gama, gsl_matrix * JJ, int N, int K, int rep, int Inum, int m, gsl_matrix * b, gsl_matrix * C,
				 char *outBetaFile, char *outZFile, char *outSigmaFile, char *outRFile, char * outGamaFile);

int main(void)
{
	int i, j, k;
	double a, rho=0.5;
	//regarding the prior of correlation R
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
	int m = 10;
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
	gsl_matrix_set_identity(Psigma);
	/*for(i = 0; i < rep; i++)
	{
		for(j = 0; j < rep; j++)
		{
			if(i == j) a = 1.0;
			else a = 0.4;
			gsl_matrix_set(Psigma, i, j, a);
		}
	}*/
	for(i=0;i<rep;i++)
	{
		for(j=0;j<rep;j++)
		{
			a = pow(rho, abs(i-j));
			gsl_matrix_set(Psigma,i,j,a);
		}
	}
	gsl_matrix_scale(Psigma, m);
	gsl_matrix_set(b, 0, 0, 1.0);
    gsl_matrix_set(b, 1, 0, 2.0);
    gsl_matrix_set_identity(C);
    gsl_matrix_scale(C, 2.0);
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
	char * outBetaFile, * outZFile, * outSigmaFile, * outRFile, * outGamaFile;
	outBetaFile = new char [80];
    outZFile = new char [80];
    outSigmaFile = new char [80];
    outRFile = new char [80];
    outGamaFile = new char [80];

	strcpy(outBetaFile, dirName);
	strcpy(outZFile, dirName);
	strcpy(outSigmaFile, dirName);
	strcpy(outRFile, dirName);
	strcpy(outGamaFile, dirName);

	strcat(outBetaFile, "CorrOrdBeta.dat");
	strcat(outZFile, "CorrOrdZ.dat");
	strcat(outSigmaFile, "CorrOrdSig.dat");
	strcat(outRFile, "CorrOrdR.dat");
    strcat(outGamaFile, "CorrOrdGama.dat");
	// Do Gibbs Sampling based on the above data.
	//Gibbsampler(Y, X, XX, Psigma, N, K, rep, Inum, m, m0, b, C, outBetaFile, outZFile, outSigmaFile);
    Gibbsampler(Y, X, YY, XX, Psigma, Gama, JJ, N, K, rep, Inum, m, b, C, outBetaFile, outZFile, outSigmaFile, outRFile, outGamaFile);
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
	delete []outSigmaFile;
	delete []outRFile;
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
	symmetric(Ainverse);

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
// generate a random matrix from inv Wishart Distribution with df = m and sigma.
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
//get correlation matrix R from W
void CorrelationM(gsl_matrix * R, gsl_matrix * Sigma)
{
	int i, j, k;
	double a, b, c;
	k = Sigma -> size1;

	for(i = 0; i < k; i++)
	{
		for(j = 0; j < k; j++)
		{
			a = gsl_matrix_get(Sigma, i, i);
			b = gsl_matrix_get(Sigma, j, j);
			c = gsl_matrix_get(Sigma, i, j);
			gsl_matrix_set(R, i, j, c / sqrt(a * b));
		}
	}
}
// sample SD , standard deviation, the sqrt root of variance from its prior
void SampleSD(gsl_matrix * SD, gsl_matrix * R, gsl_matrix * Psigma, int m)
{
    int i, j, k;
	k = Psigma -> size1;
	double a, b, c, d;

    gsl_matrix * P = gsl_matrix_alloc(k, k);
    InvWishart(P, m + k + 1, Psigma);

    for(i = 0; i < k; i++)
        {
            for (j = 0; j < k; j++)
            {
                a = gsl_matrix_get(P, i, i);
                //printf("step in Gibbs Sampling 1: %g\n", a);
                b = gsl_matrix_get(P, j, j);
                //printf("step in Gibbs Sampling 2: %g\n", b);
                c = gsl_matrix_get(R, i, j);
                //printf("step in Gibbs Sampling 3: %g\n", c);
                d = sqrt(a) * sqrt(b) * c;
                //printf("step in Gibbs Sampling 4: %g\n", d);
                gsl_matrix_set(SD, i, j, d);
            }
        }
	free(P);
}
// sample covariance and correlation matrix
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, int m, int N)
{
    double a, b, c;
    int i, j, k, dfinv;
	k = Sigma -> size1;

    gsl_matrix ** WXB, ** T;
    WXB = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
    T = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
    gsl_matrix * Txx = gsl_matrix_calloc(k, k);
    for(i = 0; i < k; i++)
    {
        for(j = 0; j < k; j++) gsl_matrix_set(Txx, i, j, 0);
    }
	for(i = 0; i < N; i++)
	{
        WXB[i] = gsl_matrix_alloc(k, 1);
        T[i] = gsl_matrix_alloc(k, k);
    }
    for(i = 0; i < N; i++)
	{
		for(j = 0; j < k; j++)
		{
			a = gsl_matrix_get(W, i * k + j, 0);
			//gsl_matrix_set(WW[i], j, 0, a);
			b = gsl_matrix_get(XB, i * k + j, 0);
			//gsl_matrix_set(XXB[i], j, 0, b);
			c = a - b;
			gsl_matrix_set(WXB[i], j, 0, c);
		}

		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, WXB[i], WXB[i], 0.0, T[i]);

		gsl_matrix_add(Txx, T[i]);
    }
    gsl_matrix_add(Txx, Psigma);

    dfinv = N + m + k + 1;
    InvWishart(Sigma, dfinv, Txx);
    CorrelationM(R, Sigma);

    gsl_matrix_free(Txx);
	for(i = 0; i < N; i++)
	{
		gsl_matrix_free(WXB[i]);
		gsl_matrix_free(T[i]);
	}
	free(WXB);
	free(T);
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

    else if(a>0) Inv(InvC, C);

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
void Gibbsampler(gsl_matrix * Y, gsl_matrix * X, gsl_matrix *YY, gsl_matrix ** XX, gsl_matrix * Psigma, gsl_matrix ** Gama, gsl_matrix * JJ, int N, int K, int rep, int Inum, int m, gsl_matrix * b, gsl_matrix * C,
				 char *outBetaFile, char *outZFile, char *outSigmaFile, char *outRFile, char * outGamaFile)
{
	int i, k;
	double a;
	gsl_matrix * mu = gsl_matrix_alloc(K, 1);
	gsl_matrix * cov = gsl_matrix_alloc(K, K);
	gsl_matrix * XB = gsl_matrix_alloc(N * rep, 1);

	// here we define Z, beta and R to replace the previous parameters;
	gsl_matrix *newZ = gsl_matrix_alloc(N * rep, 1);
	gsl_matrix *newSig = gsl_matrix_alloc(rep, rep);
	gsl_matrix *newBeta = gsl_matrix_alloc(K, 1);
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
	for(i = 1; i <= Inum; i++)
	{
       if((i % 100) == 0)
       {
	     printf("step in Gibbs Sampling : %d\n", i);
       }
       cholbeta(mu, cov, newZ, X, XX, newSig, N, b, C);
       multinorm(newBeta, mu, cov);
       gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, newBeta, 0.0, XB);
       SampleSig(newSig, newR, newZ, XB, Psigma, m, N);
       SampleZ(newZ, newSig, XB, Y, N, newGama, JJ);
       SampleGama(N, rep, YY, newZ, newGama, JJ);
       // this part is for the output;
       // this part is for the output;
       gsl_matrix_fprintf(betafp, newBeta, "%g");
       gsl_matrix_fprintf(sigfp, newSig, "%g");
       gsl_matrix_fprintf(rfp, newR, "%g");
       for(k=0; k<rep; k++)
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
	for(k=0;k<rep;k++)
	{
		gsl_matrix_free(newGama[k]);
	}
	free(newGama);
}
