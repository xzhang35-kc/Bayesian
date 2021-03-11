/////////////////////////////////////
// Name: MCMC.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions for MCMC sampling


void cholbeta(gsl_matrix * mu, gsl_matrix * cov, gsl_matrix * Z, gsl_matrix * X,
	gsl_matrix ** XX, gsl_matrix * Sigma, const int & N, gsl_matrix * b, gsl_matrix * C);
void SampleZ(gsl_matrix * Z, gsl_matrix * sigma, gsl_matrix * XB, gsl_matrix * Y,
	const int & N, gsl_matrix ** Gama, gsl_matrix * JJ);
void SampleGama(const int & N, const int & K, const int & bc, gsl_matrix * YY, gsl_matrix * Z,
	gsl_matrix ** Gama, gsl_matrix * JJ);
/*
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma,
               const int & m, const int & N);
*/
