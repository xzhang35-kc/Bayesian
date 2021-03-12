/////////////////////////////////////
// Name: Parameters.h
// Developers: Xiao & Kui Zhang at MTU
// Copyright: 2005 -
// Description: handle parameters


void printPara(const int & N, const int & P, const int & K, const int & S, const int & m, const int & m0, const int &PD,
	const int & SMH, const int & DFMH, const long int & RS, char * dirName, char * yFile, char * xFile, char * outFilePrefix);
void setPara(const int & argc, char * argv[], int & N, int & P, int & K, int & S, int & m, int & m0, int & PD,
    int & SMH, int & DFMH, long int & RS, char * dirName, char * yFile, char * xFile, char * outFilePrefix);
