#pragma once

#define USE_MKL

#ifdef USE_MKL
void _dgetrf(int m, int n, double* a, int lda, int* ipiv);
void _dgetrs(char trans, int n, int nrhs, const double* a, int lda, const int* ipiv, double* b, int ldb);
double _dnrm2(const int n, const double* X, const int incX);
void _dgemv(char TransA, const int m, const int n,
	const double alpha, const double* A, const int lda,
	const double* X, const int incX, const double beta, double* Y, const int incY);
#else
#include "cpp_lapack.h"
void _dgetrf(int m, int n, double* a, int lda, int* ipiv)
{
	fem::common cmn;
	int info;
	cpp_lapack::dgetrf(cmn, m, n, *a, lda, *ipiv, info);
}

double _dnrm2(const int n, const double* X, const int incX)
{
	return cpp_lapack::dnrm2(n, *X, incX);
}

void _dgemv(char TransA, const int m, const int n,
	const double alpha, const double* A, const int lda,
	const double* X, const int incX, const double beta, double* Y, const int incY)
{
	fem::common cmn;
	cpp_lapack::dgemv(cmn, &TransA, m, n, alpha, *A, lda, *X, incX, beta, *Y, incY);
}

void _dgetrs(char trans, int n, int nrhs, const double* a, int lda, const int* ipiv, double* b, int ldb)
{
	fem::common cmn;
	int info;
	cpp_lapack::dgetrs(cmn, &trans, n, nrhs, *a, lda, *ipiv, *b, ldb, info);
}

#endif
