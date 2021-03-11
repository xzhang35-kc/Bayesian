//////////////////////////////////////
// Name: DataInput.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: data input, set prior


void DataInput(char *dirName, char * yFile, char * xFile, gsl_matrix * Y, gsl_matrix * YY, gsl_matrix * X,
	gsl_matrix ** XX, const int & N, const int & K, const int & P);
void FindNumCut(int & ordcut, gsl_matrix * y);
void setStrSigma(gsl_matrix * Psigma, gsl_matrix * b, gsl_matrix * C, const int & m, const int & m0, const int &PD, const int & mid);
