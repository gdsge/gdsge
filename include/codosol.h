/*
CoDoSol
	  Bellavia S., Macconi M., Pieraccini S.,
	  Constrained dogleg methods for nonlinear systems with simple bounds.
	  COMPUTATIONAL OPTIMIZATION AND APPLICATIONS, vol. 53 n. 3 (2012), pp. 771-794 - ISSN 0926-6003

	  C++ version: Wenlan Luo (luowenlan@gmail.com)
	  8/23/2019
*/

#include <cmath>
#include "string.h"

#ifdef USE_SPARSE_JACOBIAN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#else
#include <Eigen/Dense>
#endif

#ifndef MAXDIM
#define MAXDIM 2
#endif

#ifndef GDSGE_SOLVE_DIM
#define GDSGE_SOLVE_DIM MAXDIM
#endif

#define INF 1e300

namespace CoDoSol
{

	inline void vec_multiply(int n, double* x, double* y, double* z)
	{
		// Calculate z = x.*y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			z[i] = x[i] * y[i];
		}

	}

	inline void vec_divide(int n, double* x, double* y, double* z)
	{
		// Calculate z = x ./ y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			z[i] = x[i] / y[i];
		}

	}
    
    inline void A_plus_x_times_yp(int n, double* A, double* x, double* y)
    {
        // Calculate A = A + x * y', column major form
        for (int i = 0; i < n; i++)
        {
            if (x[i]!=0) {
                #pragma simd
                for (int j = 0; j < n; j++)
                {
                    A[i + j*n] += x[i] * y[j];
                }
            }
        }
    }

	inline void vec_minus(int n, double* x, double* y, double* z)
	{
		// Calculate z = x - y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			z[i] = x[i] - y[i];
		}

	}
    
    inline void vec_minus_inplace(int n, double* x, double* y)
	{
		// Calculate x = x-y;
#pragma simd
		for (int i = 0; i < n; i++)
		{
			x[i] -= y[i];
		}
	}

	inline void vec_add(int n, double* x, double* y, double* z)
	{
		// Calculate z = x + y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			z[i] = x[i] + y[i];
		}

	}

	inline void dscal(int n, double* x, double alpha, double* y)
	{
		// y = alpha*x
#pragma simd
		for (int i = 0; i < n; i++)
		{
			y[i] = alpha*x[i];
		}
	}
    
    inline void dscal_inplace(int n, double* x, double alpha)
	{
		// y = alpha*x
#pragma simd
		for (int i = 0; i < n; i++)
		{
			x[i] *= alpha;
		}
	}

	inline void daxpby(int n, double* x, double a, double* y, double b, double* z)
	{
		// z = a*x + b*y
#pragma simd
		for (int i = 0; i < n; i++)
		{
			z[i] = a*x[i] + b*y[i];
		}
	}

	inline double ddot(int n, double* x, double* y)
	{
		// a = x'*y;
		double a = 0;
#pragma simd
		for (int i = 0; i < n; i++)
		{
			a += x[i] * y[i];
		}

		return a;
	}


	inline double pow2(double x)
	{
		return x*x;
	}

	inline double sign(double x)
	{
		if (x > 0)
			return 1;
		else if (x == 0)
			return 0;
		else
			return -1;
	}


	inline void dmatrici(int n, double* x, double* grad, double* l, double* u,
		double* d, double* dsqr, double* dmsqr)
	{
		// Coleman-Li Scaling Matrix
		for (int i = 0; i < n; i++)
		{
			if (grad[i] < 0)
			{
				double diff = u[i] - x[i];
				double sqdiff = sqrt(diff);
				dmsqr[i] = 1 / sqdiff;
				dsqr[i] = sqdiff;
				d[i] = diff;
			}
			else
			{
				double diff = x[i] - l[i];
				double sqdiff = sqrt(diff);
				dmsqr[i] = 1 / sqdiff;
				dsqr[i] = sqdiff;
				d[i] = diff;
			}
		}
	}

	inline double norm(int n, double* x)
	{
		//return _dnrm2(n, x, 1);
    double s = 0;
    for (int i = 0; i < n; i++)
    {
      s += pow2(x[i]);
    }
    return sqrt(s);
	}
    
    inline double norm_square(int n, double* x)
    {
    double s = 0;
    for (int i = 0; i < n; i++)
    {
      s += pow2(x[i]);
    }
    return s;
    }

	inline double vec_min(int n, double* x)
	{
		double minVal = INF;
		for (int i = 0; i < n; i++)
		{
			if (x[i] < minVal)
				minVal = x[i];
		}
		return minVal;
	}

	inline
		double max(double x, double y)
	{
		return (x > y) ? x : y;
	}

	inline
		double min(double x, double y)
	{
		return (x < y) ? x : y;
	}

	template <typename Lambda>
	double solve(int n, double* x, double* l, double* u, int flagBroyden, double atol, double rtol, int maxit, Lambda F, double* opt_info)
	{
		int ieer = 0;
		int nridut = 0;
		int itc = 0;
		int nvf = 0;

		double fx[GDSGE_SOLVE_DIM];

		F(x, fx, 0);
		nvf++;
		double fnrm = norm(n, fx);

		// Default values form params
		int maxnf = 1000;

		const double EPS = 2.2204e-16;

		double epsilon = 100 * EPS;
		double Deltamin = sqrt(EPS);
		double stoptol = atol + rtol*fnrm;

		// INTERNAL PARAMETERS
		double r = 0.1;
		double t = 0.25;
		double w = 0.75;
		double delta1 = 0.25;
		double delta2 = 2;
		double thetal = 0.99995;
		double sigma = 0.99995;

		// iteration
		double lambda = 0;
        
        #if !defined(__WIN32__) || GDSGE_SOLVE_DIM>100
	double* grad = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* grad_old = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));

	double* d = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* dsqr = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* dmsqr = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));

	double* Gp = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* dsqrgrad = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* dgrad = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* Gdgrad = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* jdgrad = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* Gmgrad = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* pcv = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	int* ipiv = (int*) malloc(sizeof(int)*(GDSGE_SOLVE_DIM));
	double* jac = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM*GDSGE_SOLVE_DIM));
	double* jacTemp = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM*GDSGE_SOLVE_DIM));
	double* sn = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* pciv = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* alp = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* aa = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* seg = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* bb = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* Gseg = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* Gpciv = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* xpp = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* fxpp = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));

	double* deltax = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* deltaf = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
	double* jacdx = (double*) malloc(sizeof(double)*(GDSGE_SOLVE_DIM));
        #else
	double grad[GDSGE_SOLVE_DIM];
	double grad_old[GDSGE_SOLVE_DIM];

	double d[GDSGE_SOLVE_DIM];
	double dsqr[GDSGE_SOLVE_DIM];
	double dmsqr[GDSGE_SOLVE_DIM];

	double Gp[GDSGE_SOLVE_DIM];
	double dsqrgrad[GDSGE_SOLVE_DIM];
	double dgrad[GDSGE_SOLVE_DIM];
	double Gdgrad[GDSGE_SOLVE_DIM];
	double jdgrad[GDSGE_SOLVE_DIM];
	double Gmgrad[GDSGE_SOLVE_DIM];
	double pcv[GDSGE_SOLVE_DIM];
	int ipiv[GDSGE_SOLVE_DIM];
	double jac[GDSGE_SOLVE_DIM*GDSGE_SOLVE_DIM];
	double jacTemp[GDSGE_SOLVE_DIM*GDSGE_SOLVE_DIM];
	double sn[GDSGE_SOLVE_DIM];
	double pciv[GDSGE_SOLVE_DIM];
	double alp[GDSGE_SOLVE_DIM];
	double aa[GDSGE_SOLVE_DIM];
	double seg[GDSGE_SOLVE_DIM];
	double bb[GDSGE_SOLVE_DIM];
	double Gseg[GDSGE_SOLVE_DIM];
	double Gpciv[GDSGE_SOLVE_DIM];
	double xpp[GDSGE_SOLVE_DIM];
	double fxpp[GDSGE_SOLVE_DIM];

	double deltax[GDSGE_SOLVE_DIM];
	double deltaf[GDSGE_SOLVE_DIM];
	double jacdx[GDSGE_SOLVE_DIM];
        #endif
        
        for (int i = 0; i < n; i++)
		{
			grad[i] = 1;
		}

		// Initial
		double Delta = 1;
		double Deltas = Delta;
		double alpha1 = INF;

		double fnrmxpp;
        
        double normdeltax = INF;
        double broydenthresh = 1e-6;

		while (fnrm > stoptol && itc < maxit && nvf < maxnf)
		{
			itc++;
			double fnrm0 = fnrm;
			// jacobian evaluation
            if (!flagBroyden || normdeltax>broydenthresh)
            {
                F(x, fx, jac);
                nvf++;
            }
            else
            {
                F(x, fx, 0);
                nvf++;
            }
            
			memcpy(grad_old, grad, sizeof(double)*n);

			#ifdef USE_SPARSE_JACOBIAN
			Eigen::Map<Eigen::MatrixXd> jac_dense(jac, n, n);
			Eigen::SparseMatrix<double> jac_sparse;
			jac_sparse = jac_dense.sparseView();
			jac_sparse.makeCompressed();
			Eigen::Map<Eigen::MatrixXd> fx_dense(fx, n, 1);
			Eigen::Map<Eigen::MatrixXd> grad_dense(grad, n, 1);
			grad_dense = jac_sparse.transpose() * fx_dense;
			#else
			Eigen::Map<Eigen::MatrixXd> jac_dense(jac, n, n);
			Eigen::Map<Eigen::MatrixXd> fx_dense(fx, n, 1);
			Eigen::Map<Eigen::MatrixXd> grad_dense(grad, n, 1);
			grad_dense = jac_dense.transpose() * fx_dense;
			#endif

			// Calculate of the scaling matrices 
			dmatrici(n, x, grad, l, u, d, dsqr, dmsqr);

			// Trust-region type assignment
			double* G = dmsqr;
			vec_multiply(n, dsqr, grad, dsqrgrad);
			vec_multiply(n, d, grad, dgrad);
			vec_multiply(n, G, dgrad, Gdgrad);
			#ifdef USE_SPARSE_JACOBIAN
			Eigen::Map<Eigen::MatrixXd> jdgrad_dense(jdgrad, n, 1);
			jdgrad_dense = jac_sparse * Eigen::Map<Eigen::MatrixXd>(dgrad, n, 1);
			#else
			Eigen::Map<Eigen::MatrixXd> jdgrad_dense(jdgrad, n, 1);
			jdgrad_dense = jac_dense * Eigen::Map<Eigen::MatrixXd>(dgrad, n, 1);
			#endif
			vec_divide(n, grad, G, Gmgrad);

			double ndsqrgrad = norm(n, dsqrgrad);
			double nGdgrad = norm(n, Gdgrad);
			double njdgrad = norm(n, jdgrad);

			// COMPUTATION OF THE MINIMIZER (PC) OF THE QUADRATIC MODEL ALONG D*GRAD 
			double vert = pow2(ndsqrgrad / njdgrad);
			dscal(n, dgrad, -vert, pcv);

			// INITIAL TRUST REGION RADIUS
			if (itc == 1)
				Delta = 1;

			// NEWTON STEP
			#ifdef USE_SPARSE_JACOBIAN
			Eigen::Map<Eigen::MatrixXd> sn_dense(sn, n, 1);
			Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> solver;
			solver.compute(jac_sparse);
			sn_dense = solver.solve(fx_dense);
			#else
			// Allocation-free dense Newton step: fixed-max-size storage (stack, <=GDSGE_SOLVE_DIM),
			// runtime size n. Avoids the per-call heap allocation that PartialPivLU<MatrixXd>
			// incurs, which serializes under OpenMP. Works for all n<=GDSGE_SOLVE_DIM (no n==GDSGE_SOLVE_DIM guard).
			typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, GDSGE_SOLVE_DIM, GDSGE_SOLVE_DIM> GDSGE_JacFixed;
			typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, GDSGE_SOLVE_DIM, 1> GDSGE_VecFixed;
			GDSGE_JacFixed jacFixed = Eigen::Map<Eigen::MatrixXd>(jac, n, n);
			GDSGE_VecFixed fxFixed  = Eigen::Map<Eigen::MatrixXd>(fx, n, 1);
			Eigen::PartialPivLU<GDSGE_JacFixed> lu(jacFixed);
			Eigen::Map<Eigen::MatrixXd>(sn, n, 1) = lu.solve(fxFixed);
			#endif

            dscal_inplace(n, sn, -1.0);

			// PROJECTED NEWTON STEP
			for (int i = 0; i < n; i++)
			{
				// Project to the boundary
				sn[i] = max(x[i] + sn[i], l[i]) - x[i];
			}
			for (int i = 0; i < n; i++)
			{
				sn[i] = min(x[i] + sn[i], u[i]) - x[i];
			}
			double nsn = norm(n, sn);
			dscal(n, sn, max(sigma, 1 - nsn), sn);

			// TRUST-REGION STRATEGY
			double rhof = 0;
			int nridu = -1;
			while (rhof<t && Delta>Deltamin)
			{
				nridu++;

				// COMPUTATION OF THE CAUCHY POINT
				if (vert*nGdgrad > Delta)
				{
					dscal(n, dgrad, -Delta / nGdgrad, pcv);
				}

				// COMPUTATION OF THE TRUNCATED CAUCHY POINT
				memcpy(pciv, pcv, sizeof(double)*n);
				double npcv = norm(n, pcv);
                alpha1 = INF;
				for (int i = 0; i < n; i++)
				{
					if (pcv[i] != 0)
					{
						alp[i] = max(
							(l[i] - x[i]) / pcv[i],
							(u[i] - x[i]) / pcv[i]
						);
					}
					else
					{
						alp[i] = INF;
					}
					if (alp[i] < alpha1)
						alpha1 = alp[i];
				}
				if (alpha1 <= 1)
				{
					dscal(n, pcv, max(thetal, 1 - npcv)*alpha1, pciv);
				}

				// PATH COMPUTATION 
				memcpy(aa, fx, sizeof(double)*n);

				#ifdef USE_SPARSE_JACOBIAN
				Eigen::Map<Eigen::MatrixXd> aa_dense(aa, n, 1);
				aa_dense += jac_sparse * Eigen::Map<Eigen::MatrixXd>(pciv, n, 1);
				#else
				Eigen::Map<Eigen::MatrixXd> aa_dense(aa, n, 1);
				aa_dense += jac_dense * Eigen::Map<Eigen::MatrixXd>(pciv, n, 1);
				#endif

				vec_minus(n, sn, pciv, seg);

				#ifdef USE_SPARSE_JACOBIAN
				Eigen::Map<Eigen::MatrixXd> bb_dense(bb, n, 1);
				bb_dense = jac_sparse * Eigen::Map<Eigen::MatrixXd>(seg, n, 1);
				#else
				Eigen::Map<Eigen::MatrixXd> bb_dense(bb, n, 1);
				bb_dense = jac_dense * Eigen::Map<Eigen::MatrixXd>(seg, n, 1);
				#endif

				double nbb = norm(n, bb);
				double gamma;
				if (nbb != 0)
				{
					double gammamin = -ddot(n, aa, bb) / pow2(nbb);
					vec_multiply(n, G, seg, Gseg);
					double nGseg = norm(n, Gseg);
					double a = pow2(nGseg);
					vec_multiply(n, G, pciv, Gpciv);
					double nGpciv = norm(n, Gpciv);
					double b = ddot(n, Gpciv, Gseg);
					double c = pow2(nGpciv) - pow2(Delta);

					double l1 = (-b + sign(-b)*sqrt(pow2(b) - a*c)) / a;
					double l2 = c / (l1*a);
					if (gammamin > 0)
					{
						double gammaend = max(l1, l2);
						gamma = min(gammamin, gammaend);
						if (gamma > 1)
						{
							for (int i = 0; i < n; i++)
							{
								if (seg[i] != 0)
								{
									alp[i] = max(
										(l[i] - (x[i] + pciv[i])) / seg[i],
										(u[i] - (x[i] + pciv[i])) / seg[i]
									);
								}
								else
								{
									alp[i] = INF;
								}
							}
							alpha1 = vec_min(n, alp);
							gamma = min(thetal*alpha1, gamma);
						}
					}
					else
					{
						double gammaend = min(l1, l2);
						gamma = max(gammamin, gammaend);
						if (gamma < 0)
						{
							for (int i = 0; i < n; i++)
							{
								if (seg[i] != 0)
								{
									alp[i] = max(
										(l[i] - (x[i] + pciv[i])) / (-seg[i]),
										(u[i] - (x[i] + pciv[i])) / (-seg[i])
									);
								}
								else
								{
									alp[i] = INF;
								}
							}
							alpha1 = vec_min(n, alp);
							gamma = max(-thetal*alpha1, gamma);
						}
					}
				} // nbb!=0
				else
				{
					gamma = 0;
				}
				double p[GDSGE_SOLVE_DIM];
				daxpby(n, pciv, 1 - gamma, sn, gamma, p);

				// ACCURACY REQUIREMENTS
				vec_add(n, x, p, xpp);
				F(xpp, fxpp, 0);
				nvf++;
				fnrmxpp = norm(n, fxpp);
				double fxpjactp[GDSGE_SOLVE_DIM];
				memcpy(fxpjactp, fx, sizeof(double)*n);
				#ifdef USE_SPARSE_JACOBIAN
				Eigen::Map<Eigen::MatrixXd> fxpjactp_dense(fxpjactp, n, 1);
				fxpjactp_dense += jac_sparse * Eigen::Map<Eigen::MatrixXd>(p, n, 1);
				#else
				Eigen::Map<Eigen::MatrixXd> fxpjactp_dense(fxpjactp, n, 1);
				fxpjactp_dense += jac_dense * Eigen::Map<Eigen::MatrixXd>(p, n, 1);
				#endif

				rhof = (fnrm - fnrmxpp) / (fnrm - norm(n, fxpjactp));
				Deltas = Delta;
				vec_multiply(n, G, p, Gp);
				Delta = min(delta1*Delta, 0.5*norm(n, Gp));

			} // while

			if (Delta <= Deltamin && rhof < t)
			{
				nridu++;
				ieer = 3;
				break;
			}
			Delta = Deltas;
            
            // UPDATING OF JACOBIAN
            if (flagBroyden)
            {
                vec_minus(n,xpp,x,deltax);
                vec_minus(n,fxpp,fx,deltaf);
                normdeltax = norm(n,deltax);
                
				Eigen::Map<Eigen::MatrixXd> jacdx_dense(jacdx, n, 1);
				jacdx_dense = jac_dense * Eigen::Map<Eigen::MatrixXd>(deltax, n, 1);
                double deltaxnormsqr = norm_square(n, deltax);
                vec_minus_inplace(n,deltaf,jacdx);
                dscal_inplace(n, deltaf, 1.0/deltaxnormsqr);
                A_plus_x_times_yp(n,jac,deltaf,deltax);
            }
			// UPDATING OF THE ITERATE
			memcpy(x, xpp, sizeof(double)*n);
			double del = 100 * EPS;
			int inddel = 0;
			for (int i = 0; i < n; i++)
			{
				if (x[i] < l[i] + del)
				{
					x[i] = l[i] + del;
					inddel = 1;
				}
				if (x[i] > u[i] - del)
				{
					x[i] = u[i] - del;
					inddel = 1;
				}
			}

			if (inddel == 0)
			{
				memcpy(fx, fxpp, sizeof(double)*n);
			}
			else {
				F(x, fx, 0);
				nvf++;
			}

			fnrm = fnrmxpp;
			double rat = fnrm / fnrm0;
			nridut += nridu;

			// STORING AND PRINTING THE ITERATION'S SUMMARY

			// UPDATING OF THE TRUST-REGION SIZE
			if (rhof > w)
			{
				Delta = max(Delta, delta2*norm(n, Gp));
			}
		}

		F(x, fx, 0);
        nvf++;
		fnrm = norm(n, fx);
        *opt_info = nvf;
        
        #if !defined(__WIN32__) || GDSGE_SOLVE_DIM>100
		free( grad );
		free( grad_old );

		free( d );
		free( dsqr );
		free( dmsqr );

		free( Gp );
		free( dsqrgrad );
		free( dgrad );
		free( Gdgrad );
		free( jdgrad );
		free( Gmgrad );
		free( pcv );
		free( ipiv );
		free( jac );
		free( jacTemp );
		free( sn );
		free( pciv );
		free( alp );
		free( aa );
		free( seg );
		free( bb );
		free( Gseg );
		free( Gpciv );
		free( xpp );
		free( fxpp );

		free( deltax );
		free( deltaf );
		free( jacdx );
		#endif
		return fnrm;
	}
}
