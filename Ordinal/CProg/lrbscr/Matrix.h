/////////////////////////////////////
// Name: Matrix.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to matrix operations, header file

void CorrelationM(gsl_matrix * R, gsl_matrix * Sigma);
void Inv(gsl_matrix * Ainverse, gsl_matrix * A);
void Mat_I_J(gsl_matrix * B, gsl_matrix * A, const int & k);
void MatI_J(gsl_matrix * B, gsl_matrix * A, const int & k);
void Mat_IJ(gsl_matrix * B, gsl_matrix * A, const int & k);
void Mat_I(gsl_matrix * B, gsl_matrix * A, const int & k);
void symmetric(gsl_matrix * A);
double logDet(gsl_matrix * A);

