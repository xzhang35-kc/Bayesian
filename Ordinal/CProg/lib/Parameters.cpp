/////////////////////////////////////
// Name: Parameters.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: handle parameters


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


void printPara(const int & N, const int & P, const int & K, const int & S, const int & m, const int & m0, const int & PD,
	const int & SMH, const int & DFMH, const long int & RS, char * dirName, char * yFile, char * xFile, char * outFilePrefix)
{
	printf("Please check carefully if the parameters are correct.\n\n");
	printf("The following parameters are used in the program:\n");
	printf("   The number of samples is: %d \n", N);
	printf("   The number of covariates is: %d \n", P);
	printf("   The number of repeated measures is: %d \n", K);
	printf("   The number of MCMC iterations is: %d \n", S);
	printf("   The d.f. for the proposed density of correlation is: %d \n", m);
	printf("   The d.f. for the prior of correlation: %d \n", m0);
	if(PD == 1)
    {
        printf("   The structure for covariance matrix is: identity matrix!\n");
    }
    else if(PD == 2)
    {
        printf("   The structure for covariance matrix is: compound symmetry !\n");
    }
	printf("   The number of iterations for inner MH (only for PX-GSM): %d \n", SMH);
	printf("   The proposed degrees of freedom for inner MH (only for PX-GSM): %d \n", DFMH);
	printf("   The seed for random number generator is: %ld \n", RS);
	printf("   The prefix for the output files is: %s \n", outFilePrefix);
	printf("   The file directory for input/output files is: %s \n", dirName);
	printf("   The file name for the responses: %s \n", yFile);
	printf("   The file name for the covariates: %s \n", xFile);
	if(PD == 1)
    {
        printf("   The results are in res-pid- sub folder!\n");
    }
    else if(PD == 2)
    {
         printf("   The results are in res-pcs- sub folder!\n");
    }
	printf("   The prefix for the output files is: %s \n", outFilePrefix);
	printf("\n");
}

void setPara(const int & argc, char * argv[], int & N, int & P, int & K, int & S, int & m, int & m0, int & PD,
	int & SMH, int & DFMH, long int & RS, char * dirName, char * yFile, char * xFile, char * outFilePrefix)
{
	int k;
	for(k = 1; k < argc; k++)
	{
		if(strcmp(argv[k], "-N") == 0 && k + 1 < argc) // total number of samples
		{
			N = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-P") == 0 && k + 1 < argc) // number of covariates
		{
			P = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-K") == 0 && k + 1 < argc) // number of repeated measures
		{
			K = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-S") == 0 && k + 1 < argc) // number of iterations
		{
			S = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-m") == 0 && k + 1 < argc) // d.f. of proposed density of correlation
		{
			m = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-m0") == 0 && k + 1 < argc) // d.f. of prior for correlation
		{
			m0 = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-PD") == 0 && k + 1 < argc) // id for prior distribution
		{
			PD = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-SMH") == 0 && k + 1 < argc) // # iterations for inner MH
		{
			SMH = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-DFMH") == 0 && k + 1 < argc) // proposed d.f. for inner MH
		{
			DFMH = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-RS") == 0 && k + 1 < argc) // seed for random number
		{
			RS = atoi(argv[k + 1]);
		}
		else if(strcmp(argv[k], "-dirName") == 0 && k + 1 < argc) // directory for input/output files
		{
			strcpy(dirName, argv[k + 1]);
		}
		else if(strcmp(argv[k], "-yFile") == 0 && k + 1 < argc) // file name for responses
		{
			strcpy(yFile, argv[k + 1]);
		}
		else if(strcmp(argv[k], "-xFile") == 0 && k + 1 < argc) // file name for covariates
		{
			strcpy(xFile, argv[k + 1]);
		}
		else if(strcmp(argv[k], "-outPre") == 0 && k + 1 < argc) // prefix for output files
		{
			strcpy(outFilePrefix, argv[k + 1]);
		}
	}

	// check and re-setup the parameters
	if(RS < 0)
    {
        time_t seconds;
        seconds = time(NULL);
        RS = seconds;
    }
	if(N <= 0 || P <= 0 || K <= 0 || S <= 0 || m <= 0 || m0 <= 0 || PD <= 0
        || PD > 2 || SMH <= 0 || DFMH <= 0)
	{
		printf("Error for input parameters! Please Check!\n");
		exit(1);
	}
	return;
}
