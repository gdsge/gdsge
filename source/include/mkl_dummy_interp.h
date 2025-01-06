/* mkl_dummy_interp.h

   A dummy header that mimics the Intel MKL 1D interpolation API, but
   internally calls user-provided cubic or linear interpolation functions.
*/

#ifndef MKL_DUMMY_INTERP_H
#define MKL_DUMMY_INTERP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>
#ifdef USE_OMP
#include <omp.h>
#endif

/* -------------------------------------------------------------------------
   1) Define needed "MKL" constants and typedefs
   ------------------------------------------------------------------------- */

/* For MKL, these are usually int or long long, etc. We'll just pick int. */
typedef int MKL_INT;

/* We mimic some MKL constants used in your code. */
#define DF_NO_HINT         0
#define DF_NO_BC           0
#define DF_BC_NOT_A_KNOT   1
#define DF_PP_NATURAL      1
#define DF_PP_DEFAULT      0
#define DF_PP_SPLINE       0
#define DF_METHOD_STD      0
#define DF_NO_IC           0

/* 
   In MKL, DFTaskPtr is a pointer to an internal structure.
   We'll mimic that with our own struct and store it as a void*.
*/
typedef void* DFTaskPtr;

/* -------------------------------------------------------------------------
   2) Declare a small internal struct to hold interpolation parameters.
   ------------------------------------------------------------------------- */
typedef struct {
    int npoints;         /* number of x-coordinates in this 1D problem */
    const double *x;     /* pointer to the x-grid (strictly increasing) */
    int vecSize;         /* how many vector functions in 'y'?          */
    const double *y;     /* pointer to the array of data, length npoints*vecSize */

    int s_order;         /* 2 => linear, 4 => cubic, etc. */
    int s_type;          /* e.g. DF_PP_NATURAL or DF_PP_DEFAULT */
    int bc_type;         /* e.g. DF_BC_NOT_A_KNOT or DF_NO_BC */

    /* Pointer to output coefficients (the user passes this in dfdEditPPSpline1D). */
    double *coeffOut;    /* user-provided array where we store polynomial coefs. */
} MyInterpTask;


/* -------------------------------------------------------------------------
   3) References to actual interpolation routines.
   ------------------------------------------------------------------------- */

/* For cubic spline not-a-knot: 
   coefs array length = ny*(npoints-1)*4, in row-major. 
   dx, a, b, c, r are user allocated, each length = ny*npoints or (npoints-1). 
*/
#include "cubic_spline_notaknot.h"

/* For linear interpolation:
   coefs array length = ny*(npoints-1)*2, in row-major
*/
#include "linear_interp_construct.h"


/* -------------------------------------------------------------------------
   5) Dummy MKL-like functions
   ------------------------------------------------------------------------- */

/* dfdNewTask1D: create a "DFTaskPtr" from the user data. */
static int dfdNewTask1D(
    DFTaskPtr *task,           /* [out] pointer to newly created "task" */
    int npoints,               /* number of points */
    const double *x,           /* the x-grid, length npoints */
    int hint1,                 /* e.g. DF_NO_HINT, not used here        */
    int vecSize,               /* how many vector components?           */
    const double *y,           /* data array, length npoints*vecSize    */
    int hint2                  /* e.g. DF_NO_HINT, not used here        */
)
{
    /* Allocate our internal struct and fill it. */
    MyInterpTask *p = (MyInterpTask*) malloc(sizeof(MyInterpTask));
    if(!p) return -1; /* out of memory */

    p->npoints = npoints;
    p->x       = x;
    p->vecSize = vecSize;
    p->y       = y;

    p->s_order = 0;   /* not known yet, set later in dfdEditPPSpline1D */
    p->s_type  = 0;
    p->bc_type = 0;
    p->coeffOut= NULL;

    *task = (DFTaskPtr)p;
    return 0; /* success */
}


/* dfdEditPPSpline1D: specify polynomial order, boundary conditions, etc. 
   This is where we store the user’s interpolation config, 
   and also the pointer to the output-coefficient array 'lastCoeff'.
*/
static int dfdEditPPSpline1D(
    DFTaskPtr task,
    int s_order, /* e.g. 4 => cubic, 2 => linear */
    int s_type,  /* e.g. DF_PP_NATURAL, DF_PP_DEFAULT */
    int bc_type, /* e.g. DF_BC_NOT_A_KNOT, DF_NO_BC */
    int li,      /* e.g. 0 => DF_NO_IC, not used here */
    int ic_type, /* e.g. DF_NO_IC       */
    int li2,     /* e.g. 0 => DF_NO_IC, not used */
    double *lastCoeff, /* user-provided array to hold coefficients */
    int hint
)
{
    MyInterpTask *p = (MyInterpTask *)task;
    if(!p) return -1;

    p->s_order  = s_order;
    p->s_type   = s_type;
    p->bc_type  = bc_type;
    p->coeffOut = lastCoeff; /* store this pointer for the next step */
    return 0; /* success */
}


/* dfdConstruct1D: this is where we actually build the polynomial
   by calling the user’s interpolation routines. 
   In MKL, 'dfdConstruct1D' finalizes the piecewise polynomial. 
*/
static int dfdConstruct1D(
    DFTaskPtr task,
    int pp_type,    /* e.g. DF_PP_SPLINE */
    int method      /* e.g. DF_METHOD_STD */
)
{
    MyInterpTask *p = (MyInterpTask *)task;
    if(!p) return -1;

    /* Decide if we do cubic-spline or linear interpolation. */
    int order = p->s_order;  /* 4 => cubic, 2 => linear. */
    int n     = p->npoints;
    int ny    = p->vecSize;

    /* We'll do a basic check. */
    if( (order != 2) && (order != 4) ) {
        return -2; 
    }

    /* We have p->coeffOut as the place to store the final coefficients. */

    /* 
       For the 1D not-a-knot spline routine, we need arrays: a, b, c, r, dx 
       each of dimension [ny*n] or [n-1]. We'll do dynamic allocation here 
    */
    if(order == 4) {
        double *a  = (double*) calloc(ny*n,     sizeof(double));
        double *b  = (double*) calloc(ny*n,     sizeof(double));
        double *c  = (double*) calloc(ny*n,     sizeof(double));
        double *rr = (double*) calloc(ny*n,     sizeof(double));
        double *dx = (double*) calloc(n-1,      sizeof(double));
        double *slopes = (double*) calloc(ny*(n-1), sizeof(double));

        if(!a || !b || !c || !rr || !dx || !slopes) {
            free(a); free(b); free(c); free(rr); free(dx); free(slopes);
            return -3; /* out of memory */
        }

        /* Call your function */
        construct_cubic_spline_notaknot(
            p->x,         /* x[]        */
            n,            /* npoints    */
            p->y,         /* y[], size n*ny */
            ny,
            p->coeffOut,  /* output coefs, length ny*(n-1)*4 */
            a, b, c, rr,  /* length ny*n each */
            dx, slopes    /* length n-1 */
        );

        free(a);
        free(b);
        free(c);
        free(rr);
        free(dx);
        free(slopes);
    }
    else if(order == 2) {
        /* linear interpolation */
        double *dx = (double*) calloc(n-1, sizeof(double));
        if(!dx) return -3; /* out of memory */

        construct_linear_interp(
            p->x,       /* x[]        */
            n,
            p->y,       /* y[], size n*ny */
            ny,
            p->coeffOut,/* output coefs, length ny*(n-1)*2 */
            dx
        );
        free(dx);
    }

    return 0; /* success */
}


/* dfDeleteTask: free the internal struct. */
static int dfDeleteTask(DFTaskPtr *task)
{
    if(!task) return 0;
    if(!(*task)) return 0;

    MyInterpTask *p = (MyInterpTask*)(*task);
    free(p);

    *task = NULL;
    return 0;
}

/*
  mkl_domatcopy: transpose (or copy) a matrix stored in row-major
                 with optional scaling by alpha.

  Parameters:
    ordering: 'R' => row-major
    trans:    'T' => transpose, 'N' => no transpose
    rows:     number of rows in the original matrix
    cols:     number of columns in the original matrix
    alpha:    scale factor to multiply each element
    src:      pointer to source matrix (row-major)
    lda:      leading dimension of src (i.e. number of columns if row-major)
    dst:      pointer to destination matrix (row-major)
    ldb:      leading dimension of dst (for row-major storage of the transposed data)

  Behavior for 'R','T':
    - We interpret the source as a rows-by-cols matrix in row-major,
      so the element at (i,j) is src[i*lda + j].
    - We write the transposed element into dst[j*ldb + i].
    - Multiply by alpha if alpha != 1.0.

  If alpha=1.0, you can skip the multiply for a slight performance gain.
*/
static void mkl_domatcopy(
    char ordering,
    char trans,
    int rows,    /* # of rows in src */
    int cols,    /* # of cols in src */
    double alpha,
    const double *src,
    int lda,
    double *dst,
    int ldb
)
{
    /* We handle only the case 'R','T' (row-major transpose). */
    if (ordering == 'R' && trans == 'T') 
    {
        /* Transpose from (rows x cols) to (cols x rows). */
        /* So the destination has 'cols' rows and 'rows' columns 
           in row-major sense, but we store it as (j,i). */

        /* If alpha == 1.0, skip the multiply for efficiency. */
        if (alpha == 1.0) {
            /* Parallel+SIMD-friendly transpose */
#ifdef USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < rows; i++) {
                /* We can optionally add #pragma omp simd or unroll in the inner loop. */
                for (int j = 0; j < cols; j++) {
                    /* dst element for transposed => (j, i) in row-major => j*ldb + i */
                    dst[(size_t)j * ldb + i] = src[(size_t)i * lda + j];
                }
            }
        }
        else {
            /* Same transpose, but multiply each element by alpha. */
#ifdef USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < rows; i++) {
                /* We can optionally add #pragma omp simd or unroll in the inner loop. */
                for (int j = 0; j < cols; j++) {
                    dst[(size_t)j * ldb + i] = alpha * src[(size_t)i * lda + j];
                }
            }
        }
    }
    else if (ordering == 'R' && trans == 'N')
    {
        /* Row-major, no transpose (copy + multiply if alpha != 1). */
        /* If needed for your code, you can implement it similarly: 
             for i in [0..rows)
               for j in [0..cols)
                 dst[i*ldb + j] = alpha * src[i*lda + j];
        */
        /* Example stub: */
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                dst[(size_t)i*ldb + j] = alpha * src[(size_t)i*lda + j];
            }
        }
    }
    else {
        /* If you need 'C','T' or 'C','N', implement here or raise an error. */
        /* For now, do nothing or raise a debug message. */
        // fprintf(stderr, "mkl_domatcopy: unhandled ordering/trans combination.\n");
    }
}

/* Other dummy functions */
void mkl_set_num_threads(int num_threads){
    // do nothing
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MKL_DUMMY_INTERP_H */
