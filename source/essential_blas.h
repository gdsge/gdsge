#pragma once

void _dgetrf(int m, int n, double* a, int lda, int* ipiv);
void _dgetrs(char trans, int n, int nrhs, const double* a, int lda, const int* ipiv, double* b, int ldb);
double _dnrm2(const int n, const double* X, const int incX);
void _dgemv(char TransA, const int m, const int n,
	const double alpha, const double* A, const int lda,
	const double* X, const int incX, const double beta, double* Y, const int incY);
