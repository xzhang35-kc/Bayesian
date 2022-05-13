/////////////////////////////////////
// Name: PX-GS-main.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: Parameter-Extended Random Walk Algorithm

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>

#include "MCMC.h"
#include "Random.h"
#include "Matrix.h"
#include "BMIDA.h"
#include "GS.h"

gsl_rng * r;


int main(int argc, char * argv[])
{
    // set up program parameters
    BMIDA bmida(argc, argv);

    // display parameters
    bmida.PrintParas();

    // allocate memory
    bmida.AllocateMemory();

    // input data sets: Y and X
    bmida.InputData();

    // set up priors
    bmida.SetPriors(1);

    // find number of cut-points
    bmida.FindNumCuts();

    // set up random number generator
   const gsl_rng_type * T;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);
   // set up random seed
   gsl_rng_set(r, bmida.RS);

   // initialize
   bmida.Initialize();

   // set output files
   bmida.SetOutFiles("GS");

   // Gibbs sampler
   bmida.GibbsGS();
}

// Gibbs sampler, main function
void BMIDA::GibbsGS()
{
   // open output files
   FILE * betafp, * zfp, * rfp, * sigfp, * gamafp;
   if((betafp = fopen(outBetaFile, "w")) == NULL)
   {
      printf("can not open file %s in GibbsMH()!\n", outBetaFile);
      exit(1);
   }
   if((zfp = fopen(outZFile, "w")) == NULL)
   {
      printf("can not open file %s in GibbsMH()!\n", outZFile);
      exit(1);
   }
   if((rfp = fopen(outRFile, "w")) == NULL)
   {
      printf("can not open file %s in GibbsMH()!\n", outRFile);
      exit(1);
   }
   if((sigfp = fopen(outSigmaFile, "w")) == NULL)
   {
      printf("can not open file %s in GibbsMH()!\n", outSigmaFile);
      exit(1);
   }
   if((gamafp = fopen(outGamaFile,"w")) == NULL)
   {
      printf("can not open file %s in GibbsMH()!\n", outGamaFile);
      exit(1);
   }

   // start Gibbs sampler
   // all matrix indexed by s - 1 are replaced by initMatrix;
   // all matrix indexed by s are replaced by newMatrix;
   for(int s = 1; s <= S; s++)
   {
      if(s % 100 ==0)
      {
          printf("steps in Gibbs sampling : %d\n", s);
      }

        // sample coefficients
      cholbeta(mu, cov, newZ, XX, newSig, N, b, InvC);
      multinorm(newBeta, mu, cov);

      // sample covariance/correlation
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, X, newBeta, 0.0, XB);
      SampleSig(newSig, newR, newZ, XB, Psigma, m0, N);

      // sample latent variables
      SampleZ(newZ, newSig, XB, Y, N, newGama, JJ, vt);

        // sample cut-points
      SampleGama(N, rep, YY, newZ, newGama, JJ, vt);

      // this part is for the output;
      gsl_matrix_fprintf(betafp, newBeta, "%g");
      gsl_matrix_fprintf(sigfp, newSig, "%g");
      gsl_matrix_fprintf(rfp, newR, "%g");
      for(int k = 0; k < rep; k++)
      {
         gsl_matrix_fprintf(gamafp, newGama[k], "%g");
      }
      if((s % 100) == 0)
            gsl_matrix_fprintf (zfp, newZ, "%g");
   }

   // close the file;
   fclose(betafp);
   fclose(zfp);
   fclose(rfp);
   fclose(sigfp);
   fclose(gamafp);
}

