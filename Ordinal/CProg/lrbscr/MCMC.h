/////////////////////////////////////
// Name: libsrc/MCMC.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to the MCMC sampling, header file

void cholbeta(gsl_matrix * mu, gsl_matrix * cov, gsl_matrix * Z, gsl_matrix ** XX,
              gsl_matrix * Sigma, const int & N, gsl_matrix * b, gsl_matrix * InvC);
void SampleZ(gsl_matrix * Z, gsl_matrix * sigma, gsl_matrix * XB, gsl_matrix * Y,
             const int & N, gsl_matrix ** Gama, gsl_matrix * JJ, const int * vt);
void SampleGama(const int & N, const int & K, gsl_matrix * YY, gsl_matrix * Z,
                gsl_matrix ** Gama, gsl_matrix * JJ, const int * vt);
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB,
               gsl_matrix * Psigma, const int & m0, const int & N);

