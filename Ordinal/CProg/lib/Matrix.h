/////////////////////////////////////
// Name: Matrix.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: functions related to matrix manipulation


void CorrelationM(gsl_matrix * R, gsl_matrix * Sigma);
void Inv(gsl_matrix * Ainverse, gsl_matrix * A);
void Mat_I_J(gsl_matrix * B, gsl_matrix * A, const int & k);
void MatI_J(gsl_matrix * B, gsl_matrix * A, const int & k);
void Mat_IJ(gsl_matrix * B, gsl_matrix * A, const int & k);
void Mat_I(gsl_matrix * B, gsl_matrix * A, const int & k);
void symmetric(gsl_matrix * A);
double Deter(gsl_matrix * A);
double Dert(gsl_matrix * A);
void Inv(gsl_matrix * Ainv, gsl_matrix * A);
void MatI_J(gsl_matrix * B, gsl_matrix * A, const int & k);
void print_gsl_matrix(gsl_matrix * M);
double print_gsl_matrix_diff(gsl_matrix * M, const int &ifp);
