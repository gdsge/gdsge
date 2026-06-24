# pragma once
/*!
 * \file InterpEval.h
 * \date 2019/08/04
 *
 * \author Wenlan Luo
 * Contact: luowenlan@gmail.com
 *
 * \brief Efficient evaluation of tensor product of splines and gradients. Support gradient evaluation
 *
 * TODO: long description
 *
 * \note
*/

#ifndef INTERP_ORDER
#define INTERP_ORDER 4
#endif

#include<string.h>

namespace InterpEval {
	template <class T>
	constexpr T pow2(T exponent) {
		return (T(1) << exponent);
	}

	inline void dscal(int n, double a, double* x)
	{
		// Standard scal: x = a*x
#pragma simd
		for (int i = 0; i < n; i++)
		{
			x[i] *= a;
		}
	}

	inline void dscal_out(int n, double a, double* x, double* y)
	{
		// Standard scal: y = a*x
#pragma simd
		for (int i = 0; i < n; i++)
		{
			y[i] = x[i]*a;
		}
	}

	inline void memcpy_double(double* dest, double* src, int size)
	{
#pragma simd
		for (int i = 0; i < size; i++)
		{
			dest[i] = src[i];
		}
	}


	inline void daxpy(int n, double a, double* x, double* y)
	{
		// Standard axpy: y = a*x + y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			y[i] = a*x[i] + y[i];
		}
	}

	inline void daxpy_out(int n, double a, double* x, double* y, double* z)
	{
		// Standard axpy: z = a*x + y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			z[i] = a*x[i] + y[i];
		}
	}

	inline void daxpby(int n, double a, double* x, double b, double* y)
	{
		// Standard axpy: y = a*x + b*y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			y[i] = a*x[i] + b*y[i];
		}
	}

	inline void daxpby_out(int n, double a, double* x, double b, double* y, double* z)
	{
		// Standard axpy: z = a*x + b*y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			z[i] = a*x[i] + b*y[i];
		}
	}

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

	template<int INTERP_EVAL_XDIM>
	inline int eval_vec(double* xSiteToLeft, int* cellOfSite, double* coefs, int fullVecEvalCoefsLength, int* order, int* pieces, double* coefsPiece)
	{
		const int xDim = INTERP_EVAL_XDIM;
		// Find pieces of coefs based on cell of site
		int offSet = 0;
		int cumProdPieces = 1;
		for (int i_dim = 0; i_dim < xDim; i_dim++)
		{
			offSet += cellOfSite[i_dim] * cumProdPieces;
			cumProdPieces *= pieces[i_dim];
		}
		// Account for order and vector length
		offSet *= fullVecEvalCoefsLength;

		// Copy coefs at offset to coefs Piece
		memcpy_double(coefsPiece, coefs + offSet, fullVecEvalCoefsLength);

		// Evaluation from the last dimension
		int evalOffset = 0;
		int increOffsetPerEval = fullVecEvalCoefsLength;
		for (int i_dim = xDim - 1; i_dim >= 0; i_dim--)
		{
			increOffsetPerEval /= INTERP_ORDER;
			// MKL stores coefficients from highest order to the lowest
			evalOffset += increOffsetPerEval*(INTERP_ORDER - 1);
			// Collapse one dimension at a time
#pragma unroll(4)
			for (int i_order = 0; i_order < INTERP_ORDER - 1; i_order++)
			{
				daxpy(increOffsetPerEval, xSiteToLeft[i_dim], coefsPiece + evalOffset, coefsPiece + evalOffset - increOffsetPerEval);
				evalOffset -= increOffsetPerEval;
			}
			// At the end of this loop, evalOffset increments by increOffsetPerEval*(order[i_dim]-1)
			// This automatically puts evalOffset pointed to the beginning for the next dim
		}

		return evalOffset;
	}

	template<int INTERP_EVAL_XDIM>
	inline int eval_with_grad_vec(double* xSiteToLeft, int* cellOfSite, double* coefs,
		int fullVecEvalCoefsLength, int* order, int* pieces, double* coefsPiece, double* grad, int* gradOffset)
	{
		const int xDim = INTERP_EVAL_XDIM;
		// Find pieces of coefs based on cell of site
		int offSet = 0;
		int cumProdPieces = 1;
		for (int i_dim = 0; i_dim < xDim; i_dim++)
		{
			offSet += cellOfSite[i_dim] * cumProdPieces;
			cumProdPieces *= pieces[i_dim];
		}
		// Account for order and vector length
		offSet *= fullVecEvalCoefsLength;

		// Copy coefs at offset to coefs Piece
		memcpy_double(coefsPiece, coefs + offSet, fullVecEvalCoefsLength);
		// Keep a pointer to data
		// double* coefs0 = coefs + offSet;
		int evalOffset = 0;

		// Set grad pointer
		double* gradAtDim[INTERP_EVAL_XDIM];
		for (int i_dim = 0; i_dim < xDim; i_dim++)
		{
			gradAtDim[i_dim] = grad + i_dim*fullVecEvalCoefsLength;
			memcpy_double(gradAtDim[i_dim], coefsPiece, fullVecEvalCoefsLength);
			gradOffset[i_dim] = 0;
		}

		// Evaluation from the last dimension
		int increOffsetPerEval = fullVecEvalCoefsLength;
		// Start from the last dimension
		for (int i_dim = xDim - 1; i_dim >= 0; i_dim--)
		{
#if INTERP_ORDER==4
			increOffsetPerEval /= 4;
			// MKL stores coefficients from highest order to the lowest
			evalOffset += increOffsetPerEval * 3;
			for (int j_dim = 0; j_dim < xDim; j_dim++)
			{
				gradOffset[j_dim] += increOffsetPerEval * 3;
			}

			// Self dimension. Multiply by the highest order
			dscal(increOffsetPerEval, 3.0, gradAtDim[i_dim] + gradOffset[i_dim]);

			// Collapse one dimension at a time, manually unroll
			// i_order = 0
			for (int i_order = 0; i_order < INTERP_ORDER - 1; i_order++)
			{
				daxpy(increOffsetPerEval, xSiteToLeft[i_dim], coefsPiece + evalOffset, coefsPiece + evalOffset - increOffsetPerEval);
				evalOffset -= increOffsetPerEval;
				// Manually loop to calcualte gradients; compiler should optimize this
				for (int j_dim = 0; j_dim < i_dim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
				{
					if (i_order < INTERP_ORDER - 2)
					{
						const int j_dim = i_dim;
						daxpby(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
							INTERP_ORDER - i_order - 2, gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
						gradOffset[j_dim] -= increOffsetPerEval;
					}
				}
				for (int j_dim = i_dim + 1; j_dim < xDim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
			}

			/*
			// i_order = 1
			{
				daxpy(increOffsetPerEval, xSiteToLeft[i_dim], coefsPiece + evalOffset, coefsPiece + evalOffset - increOffsetPerEval);
				evalOffset -= increOffsetPerEval;
				// Manually loop to calcualte gradients; compiler should optimize this
				for (int j_dim = 0; j_dim < i_dim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
				{
					const int j_dim = i_dim;
					daxpby(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						1.0, gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
				for (int j_dim = i_dim + 1; j_dim < xDim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
			}

			// i_order = 2
			{
				daxpy(increOffsetPerEval, xSiteToLeft[i_dim], coefsPiece + evalOffset, coefsPiece + evalOffset - increOffsetPerEval);
				evalOffset -= increOffsetPerEval;

				// Manually loop to calcualte gradients; compiler should optimize this
				for (int j_dim = 0; j_dim < i_dim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
				for (int j_dim = i_dim + 1; j_dim < xDim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
			}
			*/
			// At the end of this loop, evalOffset decrements by increOffsetPerEval*(order[i_dim]-1)
			// This automatically puts evalOffset pointed to the beginning for the next dim
#elif INTERP_ORDER==2
			// i_order = 0
			increOffsetPerEval /= 2;
			// MKL stores coefficients from highest order to the lowest
			evalOffset += increOffsetPerEval;
			for (int j_dim = 0; j_dim < xDim; j_dim++)
			{
				gradOffset[j_dim] += increOffsetPerEval;
			}

			// Self dimension. Multiply by the highest order
			// dscal(increOffsetPerEval, 1.0, gradAtDim[i_dim] + gradOffset[i_dim]);
			{
				daxpy(increOffsetPerEval, xSiteToLeft[i_dim], coefsPiece + evalOffset, coefsPiece + evalOffset - increOffsetPerEval);
				evalOffset -= increOffsetPerEval;
				// Manually loop to calcualte gradients; compiler should optimize this
				for (int j_dim = 0; j_dim < i_dim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
				for (int j_dim = i_dim + 1; j_dim < xDim; j_dim++)
				{
					// Other dimension. Multiply by 1  and add to grad
					daxpy(increOffsetPerEval, xSiteToLeft[i_dim], gradAtDim[j_dim] + gradOffset[j_dim],
						gradAtDim[j_dim] + gradOffset[j_dim] - increOffsetPerEval);
					gradOffset[j_dim] -= increOffsetPerEval;
				}
			}
#endif
		}

		return evalOffset;
		}


		template<int INTERP_EVAL_XDIM, int INTERP_EVAL_NUMVEC>
		inline void search_eval_vec(int numVec, double* xSite, int xDim, double** xGrid, int* xPts,
			double* coefs, int fullVecEvalCoefsLength, int* order, int* pieces, double* result, int* cellOfSite)
		{
			double xSiteToLeft[INTERP_EVAL_XDIM];
			double evalSpace[pow2<int>(2 * INTERP_EVAL_XDIM)*INTERP_EVAL_NUMVEC];
			for (int j = 0; j < xDim; ++j)
			{
				hunt(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
				xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
			}
			int evalOffset = eval_vec<INTERP_EVAL_XDIM>(xSiteToLeft, cellOfSite, coefs, fullVecEvalCoefsLength, order, pieces, evalSpace);

			memcpy_double(result, evalSpace + evalOffset, numVec);
		}

		// This avoids extra memcpy
		template<int INTERP_EVAL_XDIM, int INTERP_EVAL_NUMVEC>
		inline int search_eval_with_grad_vec_with_space(int numVec, double* xSite, int xDim, double** xGrid, int* xPts,
			double* coefs, int fullVecEvalCoefsLength, int* order, int* pieces, double* evalSpace, double* gradSpace, int* gradOffset)
		{
			double xSiteToLeft[INTERP_EVAL_XDIM];
			int cellOfSite[INTERP_EVAL_XDIM];
			for (int j = 0; j < xDim; ++j)
			{
				locate(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
				xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
			}
			int evalOffset = eval_with_grad_vec<INTERP_EVAL_XDIM>(xSiteToLeft, cellOfSite, coefs, fullVecEvalCoefsLength, order, pieces, evalSpace, gradSpace, gradOffset);
			return evalOffset;
		}

		template<int INTERP_EVAL_XDIM, int INTERP_EVAL_NUMVEC>
		inline void search_eval_with_grad_vec(double* xSite, double** xGrid, int* xPts,
			double* coefs, int fullVecEvalCoefsLength, int* order, int* pieces, double* result, double* grad, int* cellOfSite)
		{
			const int xDim = INTERP_EVAL_XDIM;
			const int numVec = INTERP_EVAL_NUMVEC;
			double xSiteToLeft[INTERP_EVAL_XDIM];
			double evalSpace[pow2<int>(2 * INTERP_EVAL_XDIM)*INTERP_EVAL_NUMVEC];
			double gradSpace[pow2<int>(2 * INTERP_EVAL_XDIM)*INTERP_EVAL_XDIM*INTERP_EVAL_NUMVEC];
			int gradOffset[INTERP_EVAL_XDIM];
			for (int j = 0; j < INTERP_EVAL_XDIM; ++j)
			{
				hunt(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
				xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
			}
			int evalOffset = eval_with_grad_vec<INTERP_EVAL_XDIM>(xSiteToLeft, cellOfSite, coefs, fullVecEvalCoefsLength, order, pieces, evalSpace, gradSpace, gradOffset);

			memcpy_double(result, evalSpace + evalOffset, numVec);
			for (int i_dim = 0; i_dim < xDim; i_dim++)
			{
				memcpy_double(grad + i_dim*numVec, gradSpace + i_dim*fullVecEvalCoefsLength + gradOffset[i_dim],
					numVec);
			}
		}
}

#undef INTERP_ORDER

