/////////////////////////////////////
// Name: lib-mh/MH.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to the MH method, header file

// void Transform(gsl_matrix * RR, gsl_matrix * WR, gsl_matrix * JJ);
void Transform(gsl_matrix * RR, gsl_matrix * WR, const int * vt);
double Priornew(const int & m0, gsl_matrix * sigma, gsl_matrix * W);
double Posterior(const int & N, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB);
int MHsampler(gsl_matrix * WR, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
   const int & N, const int & m, const int & m0, gsl_matrix * JJ, const int * vt);

