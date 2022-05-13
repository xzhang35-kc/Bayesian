/////////////////////////////////////
// Name: libsrc/BMIDA.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: BMIDA class

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <dir.h>
#include <gsl/gsl_randist.h>

extern gsl_rng * r;

#include "BMIDA.h"
#include "Random.h"

// constructor
BMIDA::BMIDA(const int & argc, char * argv[])
{
   // memory allocation for files
   dirName = new char [200];
   strcpy(dirName, "");
   xFile = new char [200];
   strcpy(xFile, "");
   yFile = new char [200];
   strcpy(yFile, "");
   outFilePrefix = new char [200];
   strcpy(outFilePrefix, "");
   outBetaFile = new char [200];
   outZFile = new char [200];
   outRFile = new char [200];
   outSigmaFile = new char [200];
   outGamaFile = new char [200];

   // initialize parameters
   N = -1;
   P = -1;
   rep = -1;
   K = -1;
   S = -1;
   m = 1500;
   m0 = 20;
   PD = 1;
   SMH = 1000;
   RS = -1;
   vt = new int [100];

   int k;
   char * vts = new char [100];
   strcpy(vts, "");
   for (k = 1; k < argc; k++)
   {
      if (strcmp(argv[k], "-N") == 0 && k + 1 < argc) // total number of samples
      {
         N = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-P") == 0 && k + 1 < argc) // number of covariates
      {
         P = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-K") == 0 && k + 1 < argc) // number of repeated measures
      {
         K = atoi(argv[k + 1]);
         rep = K;
      }
      else if (strcmp(argv[k], "-S") == 0 && k + 1 < argc) // number of iterations
      {
         S = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-m") == 0 && k + 1 < argc) // d.f. of proposed density of correlation
      {
         m = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-m0") == 0 && k + 1 < argc) // d.f. of prior for correlation
      {
         m0 = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-PD") == 0 && k + 1 < argc) // id for prior distribution
      {
         PD = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-SMH") == 0 && k + 1 < argc) // # iterations for inner MH
      {
         SMH = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-RS") == 0 && k + 1 < argc) // seed for random number
      {
         RS = atoi(argv[k + 1]);
      }
      else if (strcmp(argv[k], "-dirName") == 0 && k + 1 < argc) // directory for input/output files
      {
         strcpy(dirName, argv[k + 1]);
      }
      else if (strcmp(argv[k], "-yFile") == 0 && k + 1 < argc) // file name for responses
      {
         strcpy(yFile, argv[k + 1]);
      }
      else if (strcmp(argv[k], "-xFile") == 0 && k + 1 < argc) // file name for covariates
      {
         strcpy(xFile, argv[k + 1]);
      }
      else if (strcmp(argv[k], "-outPre") == 0 && k + 1 < argc) // prefix for output files
      {
         strcpy(outFilePrefix, argv[k + 1]);
      }
      else if (strcmp(argv[k], "-vt") == 0 && k + 1 < argc) // type of variable
      {
          strcpy(vts, argv[k + 1]);
      }
   }

   // check and re-setup the parameters
   if (RS < 0)
   {
      time_t seconds;
      seconds = time(NULL);
      RS = seconds;
   }
   if (N <= 0 || P <= 0 || K <= 0 || S <= 0 || m <= 0 || m0 <= 0 || PD <= 0
        || PD > 2 || SMH <= 0)
   {
      printf("Error for input parameters! Please Check!\n");
      exit(1);
   }
   if (strlen(xFile) == 0 || strlen(yFile) == 0)
   {
      printf("Data files for responses and covariates are needed!\n");
      exit(1);
   }

   // assign type of variables
   if (int (strlen(vts)) == 0)
   {
         for (int k = 0; k < 100; k++)
         {
            vt[k] = 0;
         }
   }
   else if (int (strlen(vts)) == 2 * K - 1)
   {
      char * tmp = new char [10];
      for (int k = 0; k < K; k++)
      {
         tmp[0] = vts[2 * k];
         tmp[1] = '\0';
         vt[k] = atoi(tmp);
         // printf("%d %s %d\n", int (strlen(tmp)), tmp, vt[k]);
         if (vt[k] != 0 && vt[k] != 1 && vt[k] != 2)
         {
            printf("variable type can only be 0, 1, or 2\n");
            printf("current argument for -vt is %s\n", vts);
            printf("please re-configue argument for -vt\n");
            exit(1);
         }
      }
      delete []tmp;
   }
   else
   {
      printf("argument %s for -vt (variable type) is incorrect!\n", vts);
      printf("please re-configue argument for -vt\n");
      exit(1);
   }

   // free memory
   delete []vts;

   return;
}

// display parameters
void BMIDA::PrintParas()
{
   printf("Please check carefully if the parameters are correct.\n\n");
   printf("The following parameters are used in the program:\n\n");
   printf("                         The number of samples is :  %d \n", N);
   printf("                      The number of covariates is :  %d \n", P);
   printf("               The number of repeated measures is :  %d \n", K);
   printf("                 The number of MCMC iterations is :  %d \n", S);
   printf("             The d.f. for the proposed covariance :  %d \n", m);
   printf("             The d.f. for the prior of covariance :  %d \n", m0);
   printf("   Number of iterations for inner MH (for PX-GSM) :  %d \n", SMH);
   printf("          The seed for random number generator is :  %ld \n", RS);
   printf("               The file name for the responses is :  %s \n", yFile);
   printf("              The file name for the covariates is :  %s \n", xFile);
   printf("                                The variables  are:  ");
   for (int k = 0; k <K; k++)
   {
      if (vt[k] == 0)
         printf("cont.");
      else if (vt[k] == 1)
        printf("ordinal");
      else
        printf("nominal");
      if (k != K - 1)
        printf(",");
      else
        printf("\n");
   }
   if (PD == 1)
   {
      printf("           The structure for covariance matrix is :  identity matrix\n");
   }
   else if (PD == 2)
   {
      printf("           The structure for covariance matrix is :  compound symmetry\n");
   }
   if (strlen(dirName) > 0)
   {
       printf("             The folder for input/output files is :  %s \n", dirName);
   }

   if (PD == 1)
   {
      printf("                    The results are in sub folder :  res-pid-\n");
   }
   else if (PD == 2)
   {
      printf("                   The results are in sub folder :  res-pcs- \n");
   }
   if (strlen(outFilePrefix) > 0)
   {
      printf("               The prefix for the output files is :  %s \n", outFilePrefix);
   }
   printf("\n");
}

void BMIDA::InputData()
{
   char * fileName = new char [200];
   FILE * fp;

   // input responses
   strcpy(fileName, dirName);
   strcat(fileName, yFile);
    if ((fp = fopen(fileName, "r")) == NULL)
   {
      printf("can not open file %s!\n", fileName);
      exit(1);
   }
   gsl_matrix_fscanf(fp, Y);
   fclose(fp);

   double a;
   for (int n = 0; n < N; n++) {
      for (int k = 0; k < K; k++) {
         a = gsl_matrix_get(Y, n * K + k, 0);
         gsl_matrix_set(YY, n, k, a);
      }
   }

   // input covariates
   strcpy(fileName, dirName);
   strcat(fileName, xFile);
    if ((fp = fopen(fileName, "r")) == NULL)
   {
      printf("can not open file %s!\n", fileName);
      exit(1);
   }
   gsl_matrix_fscanf(fp, X);
   fclose(fp);

   for (int n = 0; n < N; n++)
   {
      for (int k = 0; k < K; k++)
      {
         for (int p = 0; p < P; p++)
         {
            a = gsl_matrix_get(X, n * K + k, p);
            gsl_matrix_set(XX[n], k, p, a);
         }
      }
   }

   delete []fileName;

   return;
}

void BMIDA::SetPriors(const int & mid)
{
    // set up the scale matrix for the prior of covariance
 	if(PD == 1)  // identity matrix
	{
        gsl_matrix_set_identity(Psigma);
	}
	else if(PD == 2)  // compound symmetry; not used;
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
	}
    if(mid == 1 || mid == 2)
        gsl_matrix_scale(Psigma, m0);
    else
        gsl_matrix_scale(Psigma, 1.0 / m0);

    // set up mean for the prior of regression coefficients
    gsl_matrix_set_zero(b);

    // set up covariance for the prior of regression coefficients
    gsl_matrix_set_zero(C);
    gsl_matrix_set_zero(InvC);
}

void BMIDA::FindNumCuts()
{
   int maxNumCuts = -1;
   int a, b;
   for (int k = 0; k < K; k++)
   {
      if(vt[k] == 0)
      {
         b = 1;
         if(maxNumCuts < 1)
         {
            maxNumCuts = 1;
         }
         gsl_matrix_set(JJ, k, 0, b);
         continue;
      }
      b = -1;
      for (int n = 0; n < N; n++)
      {
         a = gsl_matrix_get(YY, n, k) + 0.000001;
         if (maxNumCuts < a + 1 && a != 999)
         {
            maxNumCuts = a + 1;
         }
         if (b < a && a != 999)
         {
            b = a;
         }
      }
      b += 1;
      if (b <= 0)
      {
         b = 1;
      }
      gsl_matrix_set(JJ, k, 0, b);
   }

   if (maxNumCuts <= 1)
   {
      printf("Check responses - all responses are continuous? \n");
      exit(1);
   }
   return;
}

void BMIDA::SetOutFiles(const char * gs)
{
   char * outDir = new char [200];
   strcpy(outDir, dirName);

   // set up output files
   if (PD == 1 && strcmp(gs, "MH") == 0)
   {
      strcat(outDir, "res-pid-px-mh/");
   }
   else if (PD == 2 && strcmp(gs, "MH") == 0)
   {
      strcat(outDir, "res-pcs-px-mh/");
   }
    else if (PD == 1 && strcmp(gs, "GS") == 0)
   {
      strcat(outDir, "res-pid-px-gs/");
   }
   else if (PD == 2 && strcmp(gs, "GS") == 0)
   {
      strcat(outDir, "res-pcs-px-gs/");
   }
   else if (PD == 1 && strcmp(gs, "GSM") == 0)
   {
      strcat(outDir, "res-pid-px-gsm/");
   }
   else if (PD == 2 && strcmp(gs, "GSM") == 0)
   {
      strcat(outDir, "res-pcs-px-gsm/");
   }
   else
   {
      strcat(outDir, "");
   }
   mkdir(outDir);

   strcpy(outBetaFile, outDir);
   strcat(outBetaFile, outFilePrefix);
   strcat(outBetaFile, "CorrOrdBeta.dat");

   strcpy(outZFile, outDir);
   strcat(outZFile, outFilePrefix);
   strcat(outZFile, "CorrOrdZ.dat");

   strcpy(outRFile, outDir);
   strcat(outRFile, outFilePrefix);
   strcat(outRFile, "CorrOrdR.dat");

   strcpy(outSigmaFile, outDir);
   strcat(outSigmaFile, outFilePrefix);
   strcat(outSigmaFile, "CorrOrdSig.dat");

   strcpy(outGamaFile, outDir);
   strcat(outGamaFile, outFilePrefix);
   strcat(outGamaFile, "CorrOrdGama.dat");

   delete []outDir;
}

void BMIDA::Initialize()
{
    // initialize covariance and correlation matrix
   gsl_matrix_set_identity(newSig);
   gsl_matrix_set_identity(newR);

   // initialize cut-points
   newGama = (gsl_matrix **) calloc(rep, sizeof(gsl_matrix *));
   for (int k = 0; k < rep; k++)
   {
      int ll = gsl_matrix_get(JJ, k, 0) + 0.000001;
      newGama[k] = gsl_matrix_calloc(ll, 1);
      if (ll > 1)
      {
          gsl_matrix_set(newGama[k], ll - 1, 0, 10000.0);
      }
   }

   // initialize latent variables Z
   /*
   double bb;
    for (int n = 0; n < N * rep; n++)
   {
      bb = gsl_ran_gaussian(r, 1.0);
      gsl_matrix_set(newZ, n, 0, bb);
   }
   */
   /*
   // For MH
   // set the initial value for Z;
   gsl_matrix_memcpy(newZ, Y);
   double a, bb;
   for (int i=0; i<N*rep; i++)
   {
      a = gsl_matrix_get(Y,i,0);
      if (a==999)
      {
         bb = gsl_ran_gaussian(r, 1.0);
         gsl_matrix_set(newZ,i,0,bb);
      }
      else if (a==0)
      {
         bb = rrtruncnorm(0, 1, 0.0);
         gsl_matrix_set(newZ, i, 0, bb);
      }
   }
   */

   /*
   // For GSM
   // set the initial value for Z;
   for (int i = 0; i < N * rep; i++)
   {
      double a = gsl_matrix_get(Y, i, 0);
      if (a > 0)
      {
         gsl_matrix_set(newZ, i, 0, rltruncnorm(0.0, 1.0, 0.0));
      }
      else if (fabs(a-999)<0.01)
      {
         gsl_matrix_set(newZ, i, 0, gsl_ran_gaussian(r, 1.0));
      }
      else
      {
         gsl_matrix_set(newZ, i, 0, rrtruncnorm(0.0, 1.0, 0.0));
      }
   }
   */

   // initialize latent variables Z
   double a, bb;
   gsl_matrix_memcpy(newZ, Y);
   for (int k = 0; k < K; k++)
   {
      int ll = gsl_matrix_get(JJ, k, 0) + 0.000001;
      for (int n = 0; n < N; n++)
      {
         a = gsl_matrix_get(Y, n * K + k, 0);
         if (fabs(a - 999) <= 1e-12)
         {
            bb = gsl_ran_gaussian(r, 1.0);
            gsl_matrix_set(newZ, n * K + k, 0, bb);
         }
         else if (fabs(a - 0) <= 1e-12)
         {
            if (ll >= 2)
            {
               bb = rrtruncnorm(0, 1, 0.0);
               gsl_matrix_set(newZ, n * K + k, 0, bb);
            }
         }
      }
   }
}

void BMIDA::AllocateMemory()
{
   // for input data
   X = gsl_matrix_alloc(N * rep, P); // covariate matrix N * K by P
   XX = (gsl_matrix **) calloc(N, sizeof(gsl_matrix *));
   for (int n = 0; n < N; n++)
   {
      XX[n] = gsl_matrix_alloc(rep, P);
   }
   Y = gsl_matrix_alloc(N * rep, 1); // response N*K
   YY = gsl_matrix_alloc(N, rep);

   // for parameter in prior distributions
   b = gsl_matrix_alloc(P, 1); // prior mean for beta
   C = gsl_matrix_calloc(P, P); // prior covariance for beta
   InvC = gsl_matrix_calloc(P, P); // inverse matrix for C, for computational purpose
   Psigma = gsl_matrix_alloc(rep, rep);

   // for cut-points
   JJ = gsl_matrix_alloc(rep, 1); // number of cut-points
   // for Gibbs sampler
   mu = gsl_matrix_alloc(P, 1); // for mean of posterior mean of coefficients
   cov = gsl_matrix_alloc(P, P); // for covariance of posterior mean of coefficients
   XB = gsl_matrix_alloc(N * rep, 1);
   newZ = gsl_matrix_alloc(N * rep, 1); // for latent variables
   newBeta = gsl_matrix_alloc(P, 1); // for regression coefficients
   newGama = (gsl_matrix **) calloc(rep, sizeof(gsl_matrix *));

   // for Gibbs sampler
   newSig = gsl_matrix_alloc(rep, rep); // for covariance matrix
   newR = gsl_matrix_alloc(rep, rep); // for correlation matrix

   // for Gibbs sampler GS and GSM
   //newSig = gsl_matrix_alloc(rep, rep);
   //newR = gsl_matrix_alloc(rep, rep);
   //SD = gsl_matrix_alloc(rep, rep);
}

BMIDA::~BMIDA()
{
   // free memory
   delete []vt;

   // free memory
   delete []dirName;
   delete []yFile;
   delete []xFile;
   delete []outBetaFile;
   delete []outZFile;
   delete []outRFile;
   delete []outSigmaFile;
   delete []outGamaFile;

   // free memory
   gsl_matrix_free(X);
   for (int n = 0; n < N; n++)
   {
      gsl_matrix_free(XX[n]);
   }
   free(XX);
   gsl_matrix_free(Y);
   gsl_matrix_free(YY);

   // free memory
   gsl_matrix_free(b);
   gsl_matrix_free(C);
   gsl_matrix_free(InvC);
   gsl_matrix_free(Psigma);

   // free memory
   gsl_matrix_free(JJ);

   // free memory
   gsl_matrix_free(mu);
   gsl_matrix_free(cov);
   gsl_matrix_free(XB);
   gsl_matrix_free(newZ);
   gsl_matrix_free(newBeta);
   for (int k = 0; k < rep; k++)
   {
      gsl_matrix_free(newGama[k]);
   }
   free(newGama);

   // free memory
   gsl_matrix_free(newSig);
   gsl_matrix_free(newR);
}
