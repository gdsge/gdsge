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
#include "essential_blas.h"

#ifdef USE_SPARSE_JACOBIAN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#endif

#ifndef MAXDIM
#define MAXDIM 2
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
                for (int j = 0; i < n; i++)
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

		double fx[MAXDIM];

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
        
        #if MAXDIM>100
	double* grad = (double*) malloc(sizeof(double)*(MAXDIM));
	double* grad_old = (double*) malloc(sizeof(double)*(MAXDIM));

	double* d = (double*) malloc(sizeof(double)*(MAXDIM));
	double* dsqr = (double*) malloc(sizeof(double)*(MAXDIM));
	double* dmsqr = (double*) malloc(sizeof(double)*(MAXDIM));

	double* Gp = (double*) malloc(sizeof(double)*(MAXDIM));
	double* dsqrgrad = (double*) malloc(sizeof(double)*(MAXDIM));
	double* dgrad = (double*) malloc(sizeof(double)*(MAXDIM));
	double* Gdgrad = (double*) malloc(sizeof(double)*(MAXDIM));
	double* jdgrad = (double*) malloc(sizeof(double)*(MAXDIM));
	double* Gmgrad = (double*) malloc(sizeof(double)*(MAXDIM));
	double* pcv = (double*) malloc(sizeof(double)*(MAXDIM));
	int* ipiv = (int*) malloc(sizeof(int)*(MAXDIM));
	double* jac = (double*) malloc(sizeof(double)*(MAXDIM*MAXDIM));
	double* jacTemp = (double*) malloc(sizeof(double)*(MAXDIM*MAXDIM));
	double* sn = (double*) malloc(sizeof(double)*(MAXDIM));
	double* pciv = (double*) malloc(sizeof(double)*(MAXDIM));
	double* alp = (double*) malloc(sizeof(double)*(MAXDIM));
	double* aa = (double*) malloc(sizeof(double)*(MAXDIM));
	double* seg = (double*) malloc(sizeof(double)*(MAXDIM));
	double* bb = (double*) malloc(sizeof(double)*(MAXDIM));
	double* Gseg = (double*) malloc(sizeof(double)*(MAXDIM));
	double* Gpciv = (double*) malloc(sizeof(double)*(MAXDIM));
	double* xpp = (double*) malloc(sizeof(double)*(MAXDIM));
	double* fxpp = (double*) malloc(sizeof(double)*(MAXDIM));

	double* deltax = (double*) malloc(sizeof(double)*(MAXDIM));
	double* deltaf = (double*) malloc(sizeof(double)*(MAXDIM));
	double* jacdx = (double*) malloc(sizeof(double)*(MAXDIM));
        #else
	double grad[MAXDIM];
	double grad_old[MAXDIM];

	double d[MAXDIM];
	double dsqr[MAXDIM];
	double dmsqr[MAXDIM];

	double Gp[MAXDIM];
	double dsqrgrad[MAXDIM];
	double dgrad[MAXDIM];
	double Gdgrad[MAXDIM];
	double jdgrad[MAXDIM];
	double Gmgrad[MAXDIM];
	double pcv[MAXDIM];
	int ipiv[MAXDIM];
	double jac[MAXDIM*MAXDIM];
	double jacTemp[MAXDIM*MAXDIM];
	double sn[MAXDIM];
	double pciv[MAXDIM];
	double alp[MAXDIM];
	double aa[MAXDIM];
	double seg[MAXDIM];
	double bb[MAXDIM];
	double Gseg[MAXDIM];
	double Gpciv[MAXDIM];
	double xpp[MAXDIM];
	double fxpp[MAXDIM];

	double deltax[MAXDIM];
	double deltaf[MAXDIM];
	double jacdx[MAXDIM];
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
			_dgemv('T',
				n, n, 1.0, jac, n, fx, 1, 0, grad, 1);
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
			_dgemv('N',
				n, n, 1.0, jac, n, dgrad, 1, 0, jdgrad, 1);
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
		 	/*
			Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(jac_dense);
			sn_dense = solver.solve(fx_dense);
			*/
			#else
			memcpy(jacTemp, jac, sizeof(double)*n*n);
			memcpy(sn, fx, sizeof(double)*n);
			_dgetrf(n, n, jacTemp, n, ipiv);
			_dgetrs('N', n, 1, jacTemp, n, ipiv, sn, n);
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
				_dgemv('N',
					n, n, 1.0, jac, n, pciv, 1, 1.0, aa, 1);
				#endif

				vec_minus(n, sn, pciv, seg);

				#ifdef USE_SPARSE_JACOBIAN
				Eigen::Map<Eigen::MatrixXd> bb_dense(bb, n, 1);
				bb_dense = jac_sparse * Eigen::Map<Eigen::MatrixXd>(seg, n, 1);
				#else
				_dgemv('N',
					n, n, 1.0, jac, n, seg, 1, 0.0, bb, 1);
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
				double p[MAXDIM];
				daxpby(n, pciv, 1 - gamma, sn, gamma, p);

				// ACCURACY REQUIREMENTS
				vec_add(n, x, p, xpp);
				F(xpp, fxpp, 0);
				nvf++;
				fnrmxpp = norm(n, fxpp);
				double fxpjactp[MAXDIM];
				memcpy(fxpjactp, fx, sizeof(double)*n);
				#ifdef USE_SPARSE_JACOBIAN
				Eigen::Map<Eigen::MatrixXd> fxpjactp_dense(fxpjactp, n, 1);
				fxpjactp_dense += jac_sparse * Eigen::Map<Eigen::MatrixXd>(p, n, 1);
				#else
				_dgemv('N',
					n, n, 1.0, jac, n, p, 1, 1.0, fxpjactp, 1);
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
                
                _dgemv('N',
					n, n, 1.0, jac, n, deltax, 1, 0.0, jacdx, 1);
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
        
        #if MAXDIM>100
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
