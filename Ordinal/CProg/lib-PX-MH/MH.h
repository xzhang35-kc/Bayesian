/////////////////////////////////////
// Name: Random.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related with MH method


void Transform(gsl_matrix * RR, gsl_matrix * WR, const int & bc);
double Priornew(const int & m0, gsl_matrix * sigma, gsl_matrix * W);
double Posterior(const int & N, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB);
int MHsampler(gsl_matrix * WR, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
            const int & N, const int & m, const int & m0, const int & bc);

