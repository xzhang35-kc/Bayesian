/////////////////////////////////////
// Name: lib-gsm/GSM.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to the GSM sampling method, header file

double Priornew(const int & m0, gsl_matrix * sigma, gsl_matrix * W);
void SampleSD(gsl_matrix * SD, gsl_matrix * R, gsl_matrix * Sigma, gsl_matrix * Psigma, const int & m, const int * vt);
int MHsampler(gsl_matrix * SD, gsl_matrix * WR, gsl_matrix * R, gsl_matrix * Psigma,
   const int & N, const int & mp, const int & m);
void AdjustZ(gsl_matrix * Z, gsl_matrix * sigma, gsl_matrix * Y, const int & N, const int * vt);
