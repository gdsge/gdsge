#pragma once

#ifndef INTERP_DIM
#define INTERP_DIM 10
#endif

#include "adept_source.h"

using adept::adouble;


inline
void locate(double xx[], int n, double x, int *j)
{
	int ju, jm, jl;

	jl = 0;
	ju = n + 1;
	while (ju - jl > 1) {
		jm = (ju + jl) >> 1;
		if (x >= xx[jm])
			jl = jm;
		else
			ju = jm;
	}
	if (x == xx[1]) *j = 1;
	else if (x == xx[n]) *j = n - 1;
	else *j = jl;
}

inline
void hunt(double xx[], int n, double x, int *jlo)
{
	int jm, jhi, inc;

	if (*jlo <= 0 || *jlo > n) {
		*jlo = 0;
		jhi = n + 1;
	}
	else {
		inc = 1;
		if (x >= xx[*jlo]) {
			if (*jlo == n) return;
			jhi = (*jlo) + 1;
			while (x >= xx[jhi]) {
				*jlo = jhi;
				inc += inc;
				jhi = (*jlo) + inc;
				if (jhi > n) {
					jhi = n + 1;
					break;
				}
			}
		}
		else {
			if (*jlo == 1) {
				*jlo = 0;
				return;
			}
			jhi = (*jlo)--;
			while (x < xx[*jlo]) {
				jhi = (*jlo);
				inc <<= 1;
				if (inc >= jhi) {
					*jlo = 0;
					break;
				}
				else *jlo = jhi - inc;
			}
		}
	}
	while (jhi - (*jlo) != 1) {
		jm = (jhi + (*jlo)) >> 1;
		if (x >= xx[jm])
			*jlo = jm;
		else
			jhi = jm;
	}
	if (x == xx[n]) *jlo = n - 1;
	if (x == xx[1]) *jlo = 1;
}

template <class DBL> inline
DBL receval4(double* coefs, int* coefsSize, int xDim, int idim, DBL* xsite, int* cellOfSite, int shift)
{
	// recursive evaluation of a piece of coefficient at certain dimension
	DBL r;
	if (idim == xDim - 1) {
		// last dimension
		double* pCoefs = coefs + shift + (*cellOfSite + 1) * 4;
		r = *(--pCoefs);
		r *= *xsite;
		r += *(--pCoefs);
		r *= *xsite;
		r += *(--pCoefs);
		r *= *xsite;
		r += *(--pCoefs);
	}
	else {
		shift += (*cellOfSite + 1) * 4;
		r = receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
		r *= *xsite;
		r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
		r *= *xsite;
		r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
		r *= *xsite;
		r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
	}
	return r;
}

template <class DBL> inline
DBL search_eval4(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize, DBL* xSite, int* cellOfSite, int vecIdx)
{
	DBL xSiteToLeft[INTERP_DIM];
	for (int j = 0; j < xDim; ++j)
	{
		hunt(xGrid[j], xPts[j] - 2, value(xSite[j]), cellOfSite + j);
		xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
	}
	return receval4(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
}

template <class DBL> inline
DBL receval2(double* coefs, int* coefsSize, int xDim, int idim, DBL* xsite, int* cellOfSite, int shift)
{
	// recursive evaluation of a piece of coefficient at certain dimension
	DBL r;
	if (idim == xDim - 1) {
		// last dimension
		double* pCoefs = coefs + shift + (*cellOfSite + 1) * 2;
		r = *(--pCoefs);
		r *= *xsite;
		r += *(--pCoefs);
	}
	else {
		shift += (*cellOfSite + 1) * 2;
		r = receval2(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
		r *= *xsite;
		r += receval2(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
	}
	return r;
}

template <class DBL> inline
DBL search_eval2(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize, DBL* xSite, int* cellOfSite, int vecIdx)
{
	DBL xSiteToLeft[INTERP_DIM];
	for (int j = 0; j < xDim; ++j)
	{
		hunt(xGrid[j], xPts[j] - 2, value(xSite[j]), cellOfSite + j);
		xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
	}
	return receval2(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
}