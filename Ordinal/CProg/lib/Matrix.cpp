/////////////////////////////////////
// Name: Matrix.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to matrix manipulation


#include <gsl/gsl_linalg.h>

#include "Matrix.h"


void print_gsl_matrix(gsl_matrix * M)
{
    int m = M -> size1;
    int n = M -> size2;
    double a;
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            a = gsl_matrix_get(M, i, j);
            printf("  %8.5e", a);
        }
        printf("\n");
    }
    printf("\n");
}

double print_gsl_matrix_diff(gsl_matrix * M, const int &ifp)
{
    int m = M -> size1;
    int n = M -> size2;
    double a, b;
    double sum = 0;
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            a = gsl_matrix_get(M, i, j);
            b = gsl_matrix_get(M, j, i);
            if(ifp > 0)
                printf("  %8.5e", a - b);
            sum = sum + fabs(a - b);
        }
        if(ifp > 0)
            printf("\n");
    }
    if(ifp > 0)
        printf("diff is: %8.5e\n", sum);
    return(sum);
}

// obtain correlation matrix R from W
void CorrelationM(gsl_matrix * R, gsl_matrix * Sigma)
{
	int K = Sigma -> size1;

	double a, b, c;
	for(int i = 0; i < K; i++)
	{
		for(int j = 0; j < K; j++)
		{
			a = gsl_matrix_get(Sigma, i, i);
			b = gsl_matrix_get(Sigma, j, j);
			c = gsl_matrix_get(Sigma, i, j);
			gsl_matrix_set(R, i, j, c / sqrt(a * b));
		}
	}
}

// obtain inverse of s general square matrix (may not be symmetric
void Inv(gsl_matrix * Ainv, gsl_matrix * A)
{
	int m = A -> size1;
	int n = A -> size2;
	if(m != n)
	{
		printf("the matrix in Inv() is not a square matrix!\n");
		exit(1);
	}
	gsl_permutation * p = gsl_permutation_alloc(m);
	gsl_matrix * Q = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(Q, A);
	int signum;
    gsl_linalg_LU_decomp(Q, p, & signum);
	gsl_linalg_LU_invert(Q, p, Ainv);

	// free of memory;
	gsl_matrix_free(Q);
	gsl_permutation_free(p);
}

// obtain the kth row from matrix A with kth element deleted.
void MatI_J(gsl_matrix * B, gsl_matrix * A, const int & k)
{
	int n = A -> size2;
	for(int i = 0; i < n - 1; i++)
	{
		if(i < k)
			gsl_matrix_set(B, 0, i, gsl_matrix_get(A, k, i));
		else
			gsl_matrix_set(B, 0, i, gsl_matrix_get(A, k, i + 1));
	}
}

// obtain the kth column from matrix A with kth element deleted.
void Mat_IJ(gsl_matrix * B, gsl_matrix * A, const int & k)
{
	int m = A -> size1;
	for(int i = 0; i < m - 1; i++)
	{
		if(i < k)
			gsl_matrix_set(B, i, 0, gsl_matrix_get(A, i, k));
		else
			gsl_matrix_set(B, i, 0, gsl_matrix_get(A, i + 1, k));
	}
}

// obtain the matrix B from matrix A with the kth row deleted.
void Mat_I(gsl_matrix * B, gsl_matrix * A, const int & k)
{
	int m = A -> size1;
	int n = A -> size2;
	for(int i = 0; i < m - 1; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if(i < k)
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i, j));
			else
				gsl_matrix_set(B, i, j, gsl_matrix_get(A, i + 1, j));
		}
	}
}

// obtain the matrix B which is without kth row and kth column of A.
void Mat_I_J(gsl_matrix * B, gsl_matrix * A, const int & k)
{
	int m = A -> size1;
	int n = A -> size2;
	for(int i = 0; i < m - 1; i++)
	{
		for(int j = 0; j < n - 1; j++)
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

// make a matrix to be symmetric
// due to numerical computation, a symmetric matrix may become non-symmetric
void symmetric(gsl_matrix * A)
{
	int m = A -> size1;
	int n = A -> size2;

	if(m != n)
	{
		printf("the input matrix in symmetric() is not a square matrix!\n");
		exit(1);
	}

	double a, b, d;
	for(int i = 0; i< m; i++)
	{
		for(int j = 0; j < i; j++)
		{
			a = gsl_matrix_get(A, j, i);
			b = gsl_matrix_get(A, i, j);
			d = (a + b) / 2;
			gsl_matrix_set(A, i, j, d);
			gsl_matrix_set(A, j, i, d);
		}
	}
}

// compute log determinant of a positive definite matrix A
// In this routine, A is not changed.
double Deter(gsl_matrix * A)
{
	int K = A -> size1;
	gsl_matrix * Q = gsl_matrix_alloc(K, K);
    gsl_permutation * p = gsl_permutation_alloc(K);

	symmetric(A);
	gsl_matrix_memcpy(Q, A);
    int signum;
    gsl_linalg_LU_decomp(Q, p, &signum);
    double sum = log(gsl_linalg_LU_det(Q, signum));

	gsl_matrix_free(Q);
	gsl_permutation_free(p);

	return(sum);
}

// computing the (1/2) log determination of a positive definite matrix A.
// A is not changed.
double Dert(gsl_matrix * A)
{
	int K = A -> size1;
	gsl_matrix * Q = gsl_matrix_alloc(K, K);

	symmetric(A);
	gsl_matrix_memcpy(Q, A);

	// choleskey decomposition
	gsl_linalg_cholesky_decomp(Q);

	// get the lower triangle matrix Sig
	double a;
	double sum = 0;
	for(int i = 0; i < K; i++)
	{
		a = gsl_matrix_get(Q, i, i);
		sum = sum + log(a);
	}

	// free memory;
	gsl_matrix_free(Q);

	return(sum);
}

