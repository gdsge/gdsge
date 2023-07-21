#define MAXGRADIENT 1000

#ifndef MAXDIM
#define MAXDIM 10
#endif

#include "adept_source.h"

inline
double value(const adept::adouble& x)
{
	return x.value();
}

inline
double value(double x)
{
	return x;
}

namespace adept {
	inline void memcpy_double(double* dest, double* src, int size)
	{
#pragma simd
		for (int i = 0; i < size; i++)
		{
			dest[i] = src[i];
		}
	}

	// Extension for jacobian evaluation with vectorization
	// The compiler can now do a better job than manual multipass
#ifdef __WIN32__
	void Stack::jacobian_forward_vec(double* jacobian_out)
	{
        const Offset numIndependent = MAXDIM;
        double* gradient;
        #if MAXDIM>100
        std::vector<double> gradient_vec(max_gradient_*numIndependent);
        gradient = gradient_vec.data();
        #else
		// double gradient_stack[MAXGRADIENT*numIndependent] = { 0.0 };
		double gradient_stack[MAXGRADIENT*numIndependent];
		std::vector<double> gradient_vec;
		if (max_gradient_ > MAXGRADIENT)
		{
			// Allocate space on heap
			gradient_vec.resize(max_gradient_*numIndependent);
			gradient = gradient_vec.data();
		}
		else
		{
			// Use stack
			gradient = gradient_stack;
		}
        #endif
		memset(gradient, 0, sizeof(double)*max_gradient_*numIndependent);

		/*
		std::vector<double> gradient_vec(max_gradient_*numIndependent, 0);
		double* gradient = gradient_vec.data();
		*/
		
#define GRADIENT(m,n) gradient[(m)*numIndependent+(n)]

		// Initialize
		for (Offset i = 0; i < n_independent(); i++) {
			GRADIENT(independent_offset_[i], i) = 1.0;
		}

		for (Offset ist = 1; ist < n_statements_; ist++) {
			const Statement& statement = statement_[ist];
			double a[numIndependent];
			memset(a, 0, sizeof(double)*numIndependent);
			// double* a = &GRADIENT(statement.offset, 0);
			for (Offset iop = statement_[ist - 1].end_plus_one;
				iop < statement.end_plus_one; iop++) {
#pragma code_align 64
#pragma vector aligned
				for (Offset i = 0; i < numIndependent; i++) {
                    a[i] += multiplier_[iop] * GRADIENT(offset_[iop], i);
				}
			}
			memcpy_double(&GRADIENT(statement.offset, 0), a, numIndependent);
		}
		// Copy the gradients corresponding to the dependent variables
		// into the Jacobian matrix
		for (Offset idep = 0; idep < n_dependent(); idep++) {
			for (Offset i = 0; i < n_independent(); i++) {
				jacobian_out[i*n_dependent() + idep]
					= GRADIENT(dependent_offset_[idep], i);
			}
		}
#undef GRADIENT
	}
#endif
}
