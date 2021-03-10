/////////////////////////////////////
// Name: Random.cpp
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related with MH method

void MHWishart(gsl_matrix * w, const int & m, gsl_matrix * sigma);
void Transform(gsl_matrix * RR, gsl_matrix * WR, const int & bc);
double Dert(gsl_matrix * A);
double Priornew(const int & m0, gsl_matrix * sigma, gsl_matrix * W);
double Posterior(const int & N, gsl_matrix * R, gsl_matrix * Z, gsl_matrix * XB);
int MHsampler(gsl_matrix * WR, gsl_matrix * Z, gsl_matrix * XB, gsl_matrix * Psigma,
            const int & N, const int & m, const int & m0, const int & bc);
void MHmultinorm(gsl_matrix * Mnorm, gsl_matrix * mu, gsl_matrix * sigma);

