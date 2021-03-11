/////////////////////////////////////
// Name: DataInput.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: data input, set prior


#include <string.h>
#include <gsl/gsl_matrix.h>


void DataInput(char *dirName, char * yFile, char * xFile, gsl_matrix * Y, gsl_matrix * YY, gsl_matrix * X,
	gsl_matrix ** XX, const int & N, const int & K, const int & P)
{
	char * fileName = new char [200];
    FILE * fp;

	// input responses
    strcpy(fileName, dirName);
    strcat(fileName, yFile);
 	if((fp = fopen(fileName, "r")) == NULL)
	{
		printf("can not open file %s!\n", fileName);
		exit(1);
	}
    gsl_matrix_fscanf(fp, Y);
    fclose(fp);

	double a;
    for(int n = 0; n < N; n++)
	{
		for(int k = 0; k < K; k++)
		{
			a = gsl_matrix_get(Y, n * K + k, 0);
			gsl_matrix_set(YY, n, k, a);
		}
	}

	// input covariates
	strcpy(fileName, dirName);
    strcat(fileName, xFile);
 	if((fp = fopen(fileName, "r")) == NULL)
	{
		printf("can not open file %s!\n", fileName);
		exit(1);
	}
    gsl_matrix_fscanf(fp, X);
    fclose(fp);

	for(int n = 0; n < N; n++)
	{
		for(int k = 0; k < K; k++)
		{
			for(int p = 0; p < P; p++)
			{
				a = gsl_matrix_get(X, n * K + k, p);
				gsl_matrix_set(XX[n], k, p, a);
			}
		}
	}

	delete []fileName;

	return;
}

// find number of cut-off points
// new to revise later
void FindNumCut(int & ordcut, gsl_matrix * y)
{
	ordcut = -1;
	int a;
	int N = y -> size1;
	int K = y -> size2;
	for(int n = 0; n < N; n++)
	{
		for(int k = 0; k < K; k++)
		{
			a = gsl_matrix_get(y, n * K + k, 0) + 0.0001;
			if(ordcut < a && a != 999)
				ordcut = a;
		}
	}
	ordcut = ordcut + 1;
	if(ordcut <= 0)
	{
		printf("Check response - !!!\n");
		exit(1);
	}
	return;
}

// set structure of covariance matrix
void setStrSigma(gsl_matrix * Psigma, gsl_matrix * b, gsl_matrix * C, const int & m, const int & m0, const int &PD, const int & mid)
{
    int rep =  Psigma -> size1;
    int P = b -> size1;
	if(PD == 1)  // identity matrix
	{
        gsl_matrix_set_identity(Psigma);
        if(mid == 1 || mid == 2)
            gsl_matrix_scale(Psigma, m);
        else
            gsl_matrix_scale(Psigma, 1.0 / m0);
        for(int p = 0; p < P; p++)
        {
            gsl_matrix_set(b, p, 0, 0);
        }
        gsl_matrix_set_zero(C);
	}
	else if(PD == 2)  // identity matrix
	{
	    double a;
        gsl_matrix_set_identity(Psigma);
        for(int i = 0; i < rep; i++)
        {
            for(int j = 0; j < rep; j++)
            {
                if(i == j)
                    a = 1.0;
                else
                    a = 0.4;
                gsl_matrix_set(Psigma, i, j, a);
            }
        }
        if(mid == 1 || mid == 2)
            gsl_matrix_scale(Psigma, m);
        else
            gsl_matrix_scale(Psigma, 1.0 / m0);
        for(int p = 0; p < P; p++)
        {
            gsl_matrix_set(b, p, 0, 3);
        }
        gsl_matrix_set_identity(C);
    }
}
