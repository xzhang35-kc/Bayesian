/////////////////////////////////////
// Name: GSM.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: routines for GSM method


double Priornew(const int & m0, gsl_matrix * sigma, gsl_matrix * W);
void SampleSD(gsl_matrix * SD, gsl_matrix * R, gsl_matrix * Psigma, const int & m);
void SampleSig(gsl_matrix * Sigma, gsl_matrix * R, gsl_matrix * W, gsl_matrix * XB, gsl_matrix * Psigma, const int & m, const int & N);
int MHsampler(gsl_matrix * SD, gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
              const int & N, const int & pm, const int & m);
