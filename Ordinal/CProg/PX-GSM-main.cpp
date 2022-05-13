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
#include "GSM.h"

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
    bmida.SetOutFiles("GSM");

    // Gibbs sampler
    bmida.GibbsGSM();
}

// Gibbs sampler, main function
void BMIDA::GibbsGSM()
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

	// for MH method
    gsl_matrix * SD = gsl_matrix_alloc(rep, rep); // for correlation matrix
    gsl_matrix * WR = gsl_matrix_alloc(rep, rep); // for covariance matrix

	// initialize
	if(PD == 1)
        gsl_matrix_set_identity(SD);

	// start Gibbs sampler
	// all matrix indexed by s - 1 are replaced by initMatrix;
	// all matrix indexed by s are replaced by newMatrix;
	int sp = 0;
    for(int s = 1; s <= S; s++)
	{
	    if((s % 100) == 0)
        {
            printf("steps in Gibbs sampling : %d\n", s);
            if(PD == 2)
                printf("accepted in MH: %d\n", sp);
        }

        // sample covariance/correlation
        if(PD == 1)
        {
            SampleSD(SD, newR, newSig, Psigma, m0, vt);
        }
        else if(PD == 2)
        {
            gsl_matrix_memcpy(SD, newSig);
            gsl_matrix_memcpy(WR, newSig);
            sp = 0;
            for(int i = 0; i < SMH; i++)
                sp = sp + MHsampler(SD, WR, newR, Psigma, N, m, m0);
        }

		// sample latent variables
		SampleZ(newZ, SD, XB, Y, N, newGama, JJ, vt);
		AdjustZ(newZ, SD, Y, N, vt);

        // sample cut-points
		SampleGama(N, rep, YY, newZ, newGama, JJ, vt);

		// sample covariance/correlation
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, newBeta, 0.0, XB);
      SampleSig(newSig, newR, newZ, XB, Psigma, m0, N);

        // sample coefficients
		cholbeta(mu, cov, newZ, XX, newSig, N, b, InvC);
		multinorm(newBeta, mu, cov);

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

	gsl_matrix_free(SD);
	gsl_matrix_free(WR);
}

