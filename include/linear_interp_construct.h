/*
    linear_interp_construct.h

    A simple header for constructing linear interpolation coefficients for
    vector-valued data. The data format is similar to the cubic spline
    code in "cubic_spline_notaknot.h".

    We assume:
      - x[0..nx-1]: strictly increasing knots
      - y[0..ny*nx-1]: row-major data for ny vector functions
          For the function with index iy, samples are y[iy*nx + 0..nx-1].
      - coefs: an output array of size (ny*(nx-1)*2), also in row-major order.
         The 2 coefficients per interval are (constantTerm, slope),
         i.e. for the i-th interval of the iy-th function:
           coefs[ (iy*(nx-1) + i)*2 + 0 ] = y[i]
           coefs[ (iy*(nx-1) + i)*2 + 1 ] = (y[i+1] - y[i]) / dx[i]
      - dx: a user-allocated array of length (nx-1), which should store
          x[i+1] - x[i]. This is either precomputed by the user or can be
          filled in by this function if desired.

    This function is allocation-free. If compiled with -DUSE_OMP and linked
    against OpenMP, it can parallelize across the ny functions.

    Example usage:

      #include "linear_interp_construct.h"

      // ...
      // Suppose you have x[] of length nx, y[] of length (ny*nx).
      // You want coefs[] of length (ny*(nx-1)*2).
      // Also have dx[] of length (nx-1).

      construct_linear_interp(
          x, nx,
          y, ny,
          coefs,
          dx
      );

      // The array 'coefs' now contains the linear interpolation coefficients.
*/

#ifndef LINEAR_INTERP_CONSTRUCT_H
#define LINEAR_INTERP_CONSTRUCT_H

#include <stdlib.h>
#ifdef USE_OMP
#include <omp.h>
#endif

static void construct_linear_interp(
    const double *x,   /* length nx */
    int nx,
    const double *y,   /* length ny*nx, row-major */
    int ny,
    double *coefs,     /* length (ny*(nx-1)*2), output in row-major */
    double *dx         /* length (nx-1), assumed allocated by caller */
)
{
    if (nx < 2) return;

    /* If not precomputed, fill dx[i] = x[i+1] - x[i]. */
    for(int i = 0; i < nx - 1; i++){
        dx[i] = x[i+1] - x[i];
    }

    /* For each vector function, build the (constant, slope) pair per interval. */
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for(int iy = 0; iy < ny; iy++)
    {
        const double *yrow = &y[iy * nx];
        double *myCoefs    = &coefs[(iy*(nx-1)) * 2];
        for(int i = 0; i < nx - 1; i++)
        {
            double slope = (yrow[i+1] - yrow[i]) / dx[i];
            myCoefs[2*i + 0] = yrow[i];   /* constant term = value at left knot */
            myCoefs[2*i + 1] = slope;     /* slope */
        }
    }
}

#endif /* LINEAR_INTERP_CONSTRUCT_H */
