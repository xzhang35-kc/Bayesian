/////////////////////////////////////
// Name: libscr/Random.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to random numbers, header file

void InvWishart(gsl_matrix * SigmaMat, const int & m, gsl_matrix * ScaleMat);
void Wishart(gsl_matrix * w, const int & m, gsl_matrix * Sigma);
double rltruncnorm(const double & mu, const double & sd, const double & L);
double rrtruncnorm(const double & mu, const double & sd, const double & L);
double truncnorm(const double & mu, const double & sd, const double & a, const double & b);
void multinorm(gsl_matrix * Mnorm, gsl_matrix * mu, gsl_matrix * sigma);

