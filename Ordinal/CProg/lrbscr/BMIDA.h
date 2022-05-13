/////////////////////////////////////
// Name: libsrc/BMIDA.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: BMIDA class, head file

# include <gsl/gsl_linalg.h>

class BMIDA {
   public:
      // parameters used in the program
      int N; // number of samples
      int P; // number of covariates
      int rep; // number of repeated measures
      int K; // same as rep, number of repeated measures
      int S; // number of iterations
      int m; // d.f. for the proposed density of correlation
      int m0; // d.f. for the prior of correlation
      int PD; // prior for Sigma, default is 1 and identity
      int SMH; // number of iteration for inner MH
      long int RS; // seed for random number; default -1: a random seed
      int * vt; // type of variable: 0 - continuous; 1 - ordinal; 2 - nominal

      // files related
      char * dirName;
      char * xFile;
      char * yFile;
      char * outFilePrefix;
      char * outBetaFile;
      char * outZFile;
      char * outRFile;
      char * outSigmaFile;
      char * outGamaFile;

      // for input data
      gsl_matrix * X;
      gsl_matrix ** XX;
      gsl_matrix * Y;
      gsl_matrix * YY;

      // for parameter in prior distributions
      gsl_matrix * b;
      gsl_matrix * C;
      gsl_matrix * InvC;
      gsl_matrix * Psigma;

      // for cut-points
      gsl_matrix * JJ;

      // for Gibbs sampler
      gsl_matrix * mu; // for mean of posterior mean of coefficients
      gsl_matrix * cov; // for covariance of posterior mean of coefficients
      gsl_matrix * XB;
      gsl_matrix * newZ; // for latent variables
      gsl_matrix * newBeta; // for regression coefficients
      gsl_matrix ** newGama; // for cut-points

      // for Gibbs sampler
      gsl_matrix * newSig; // for covariance matrix, same as newSig
      gsl_matrix * newR; // for correlation matrix, same as newR

      // some functions
      void PrintParas(); // display parameters
      void AllocateMemory();
      void InputData();
      void FindNumCuts();
      void SetPriors(const int & mid);
      void SetOutFiles(const char * gs);
      void Initialize();

      // Gibbs sampler
      void GibbsMH();
      void GibbsGS();
      void GibbsGSM();

      // constructor and destructor
      BMIDA(const int & argc, char * argv[]); // constructor
      BMIDA(); // constructor
      ~BMIDA(); // destructor
};
