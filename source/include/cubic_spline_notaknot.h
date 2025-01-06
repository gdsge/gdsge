/*
    cubic_spline_notaknot.h

    A standalone C header implementing cubic spline interpolation with
    "not-a-knot" boundary conditions for vector-valued data. 

    This version allows the caller to pass in pre-allocated arrays 
    a[], b[], c[], r[] of size (ny * nx), so that each vector function
    (indexed by iy in [0..ny-1]) will use a slice of length nx in these
    arrays. This avoids re-allocations inside the function and helps 
    prevent data races when running under OpenMP.

    The calling usage is as follows:

      // The user will allocate:
      double *a = (double*) malloc(ny * nx * sizeof(double));
      double *b = (double*) malloc(ny * nx * sizeof(double));
      double *c = (double*) malloc(ny * nx * sizeof(double));
      double *r = (double*) malloc(ny * nx * sizeof(double));
      double *intervals = (double*) malloc((nx-1) * sizeof(double));
      double *slopes = (double*) malloc(ny * (nx-1) * sizeof(double));

      // Then call the spline constructor:
      construct_cubic_spline_notaknot(
          x, nx,
          y, ny,
          coefs,
          a, b, c, r,
          intervals,
          slopes
      );

      // After the call, coefs[ ny*(nx-1)*4 ] contains the
      // desired spline coefficients.

    Notation:
      - x[0..nx-1] : sorted knot locations (must be strictly increasing).
      - y[0..(ny*nx)-1] : row-major data for ny vector functions.
          For each vector function iy, the data is y[iy*nx + i].
      - coefs[0..(ny*(nx-1)*4)-1] : output array of polynomial coeffs:
          coefs[ (iy*(nx-1) + iInterval)*4 + iCoef ] in ascending 
          powers: c0, c1, c2, c3.
      - intervals[0..nx-2] : pre-allocated array for storing knot intervals
      - slopes[0..(ny*(nx-1))-1] : pre-allocated array for storing slopes

    The code internally:
      1) For each vector function, it forms a tridiagonal system
         of size nx enforcing not-a-knot boundary conditions.
      2) Solves for the "second derivative" M[i] at each knot i.
      3) Converts M[] to local polynomial coefficients on each subinterval.

    We use a tridiagonal solver with partial pivoting for improved numerical 
    stability when dealing with small diagonal elements.

    If compiled with -fopenmp, the loops over the ny vector functions
    can be parallelized automatically.

    (C) 2024, Provided as-is, no warranty.
*/

#ifndef CUBIC_SPLINE_NOTAKNOT_H
#define CUBIC_SPLINE_NOTAKNOT_H

#include <stdlib.h>
#include <math.h>    /* for fabs */
#ifdef USE_OMP
#include <omp.h>
#endif

inline void solve_tridiagonal_system(
    double* lower,
    double* diag,
    double* upper,
    double* rhs,
    int dim
)
{
    // Elimination step.
    for (int i = 0; i < dim - 1; i++) {
        if (fabs(diag[i]) >= fabs(lower[i + 1])) {
            // Normal tridiagonal algorithm.
            auto const w = lower[i + 1] / diag[i];
            diag[i + 1] -= w * upper[i];
            rhs[i + 1] -= w * rhs[i];
            lower[i + 1] = 0;
        } else {
            // Pivoting. Here, we interchange the i-th row and the (i+1)-th
            // row, then eliminate. Unlike the other branch, the lower
            // triangular element lower[i+1] will remain. This affects the
            // back substitution step below.
            auto const w = diag[i] / lower[i + 1];

            auto const u = diag[i + 1];
            diag[i] = lower[i + 1];
            diag[i + 1] = upper[i] - w * u;
            lower[i + 1] = upper[i + 1];
            upper[i + 1] *= -w;
            upper[i] = u;

            auto const r = rhs[i];
            rhs[i] = rhs[i + 1];
            rhs[i + 1] = r - w * rhs[i + 1];
        }
    }

    // Back-substitution step.
    rhs[dim - 1] /= diag[dim - 1];

    for (int i_rev = 2; i_rev <= dim; i_rev++) {
        auto const i = dim - i_rev;
        if (i == dim - 2) {
            rhs[i] -= upper[i] * rhs[i + 1];
        } else {
            rhs[i] -= upper[i] * rhs[i + 1] + lower[i + 1] * rhs[i + 2];
        }
        rhs[i] /= diag[i];
    }
}

void construct_cubic_spline_notaknot(
    const double *knots,    /* length nx */
    int nx,
    const double *y,    /* length ny*nx */
    int ny,
    double *_coefs,      /* length ny*(nx-1)*4, output */
    double *_L,          /* length ny*nx (sub diag) */
    double *_D,          /* length ny*nx (diag) */
    double *_U,          /* length ny*nx (sup diag) */
    double *_Y,          /* length ny*nx (rhs / solution) */
    double *intervals,  /* length nx-1 */
    double *_slopes      /* length ny*(nx-1) */
)
{
    auto const n_segments = nx - 1;
    for (int i = 0; i < n_segments; i++) {
      intervals[i] = knots[i + 1] - knots[i];
    }
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iy = 0; iy < ny; iy++) {
        auto const values = y + (size_t)iy * nx;
        auto const slopes = _slopes + (size_t)iy * n_segments;
        auto const L = _L + (size_t)iy * nx;
        auto const D = _D + (size_t)iy * nx;
        auto const U = _U + (size_t)iy * nx;
        auto const Y = _Y + (size_t)iy * nx;
        auto const coefs = _coefs + (size_t)iy * (nx - 1) * 4;

        for (int i = 0; i < n_segments; i++) {
            auto const dx = values[i + 1] - values[i];
            slopes[i] = dx / intervals[i];
        }

        {
            auto const h0 = intervals[0];
            auto const h1 = intervals[1];
            auto const s0 = slopes[0];
            auto const s1 = slopes[1];
            D[0] = h0 - h1;
            U[0] = 2 * h0 + h1;
            Y[0] = -6 * h0 / (h0 + h1) * (s0 - s1);
        }

        {
            auto const h0 = intervals[n_segments - 1];
            auto const h1 = intervals[n_segments - 2];
            auto const s0 = slopes[n_segments - 1];
            auto const s1 = slopes[n_segments - 2];
            D[n_segments] = h0 - h1;
            L[n_segments] = 2 * h0 + h1;
            Y[n_segments] = 6 * h0 / (h0 + h1) * (s0 - s1);
        }

        // The remaining coefficients are derived from the spline conditions.
        for (int i = 1; i < n_segments; i++) {
            L[i] = intervals[i - 1];
            D[i] = 2 * (intervals[i - 1] + intervals[i]);
            U[i] = intervals[i];
            Y[i] = 6 * (slopes[i] - slopes[i - 1]);
        }

        // Solve the equations. Solution is returned to Y.
        solve_tridiagonal_system(L, D, U, Y, n_segments + 1);
        auto const& M = Y;

        // Derive the polynomial coefficients of each spline segment from the
        // second derivatives `M`.
        for (int i = 0; i < n_segments; i++) {
            int idx = i*4;
            coefs[idx + 0] = values[i];
            coefs[idx + 1] = slopes[i] - (M[i + 1] + 2 * M[i]) * intervals[i] / 6;
            coefs[idx + 2] = M[i] / 2;
            coefs[idx + 3] = (M[i + 1] - M[i]) / (6 * intervals[i]);
        }
    }
}

#endif /* CUBIC_SPLINE_NOTAKNOT_H */
