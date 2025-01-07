/**
	Adaptive Sparse Grid Interpolation that provides adouble evaluation
	@author: Wenlan Luo, luowenlan@gmail.com
*/

#ifndef ASG_MAX_DIM
#define ASG_MAX_DIM 10
#endif 

#ifndef ASG_MAX_NVEC
#define ASG_MAX_NVEC 100
#endif

#ifndef ASG_MAX_LEVEL
#define ASG_MAX_LEVEL 20
#endif

#include "asg.h"
#include "adept_source.h"

using adept::adouble;

typedef std::array<adouble, ASG_MAX_DIM> inarray_adouble;

namespace AdaptiveSparseGridInterp {
	inline
	void search_adept(adouble* site, double* cell, adouble* ratio, int searchLevel, int nDim)
	{
#define CELL(m,n) cell[(m)*nDim+(n)]
#define RATIO(m,n) ratio[(m)*nDim+(n)]

		double pow_2_next_i_level = 1.0;
		for (int i_level = 0; i_level <= searchLevel; i_level++)
		{
			pow_2_next_i_level *= 0.5;
			if (i_level == 0)
			{
				// Special case for level 0
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					CELL(i_level, i_dim) = 0.5;
					RATIO(i_level, i_dim) = 1.0;
				}
			}
			else if (i_level == 1)
			{
				// Special case for level 1
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					site[i_dim] *= 2;
					if (site[i_dim] < 1)
					{
						CELL(i_level, i_dim) = 0.0;
						RATIO(i_level, i_dim) = 1 - site[i_dim];
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = 0.25;
					}
					else
					{
						CELL(i_level, i_dim) = 1.0;
						site[i_dim] -= 1;
						RATIO(i_level, i_dim) = site[i_dim];
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = 0.75;
					}
				} // i_dim
			} // i_level==1
			else
			{
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					site[i_dim] *= 2;
					if (site[i_dim] < 1)
					{
						RATIO(i_level, i_dim) = site[i_dim];
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = CELL(i_level, i_dim) - pow_2_next_i_level;
					}
					else
					{
						site[i_dim] -= 1;
						RATIO(i_level, i_dim) = (1 - site[i_dim]);
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = CELL(i_level, i_dim) + pow_2_next_i_level;
					}
				}
			} // i_level != 1
		}
	}

	inline
	void search_adept_with_slope(double* site, double* cell, double* ratio, double* slope, int searchLevel, int nDim)
	{
#define CELL(m,n) cell[(m)*nDim+(n)]
#define RATIO(m,n) ratio[(m)*nDim+(n)]
#define SLOPE(m,n) slope[(m)*nDim+(n)]

		double pow_2_next_i_level = 1.0;
		double twoToPowerLevel = 1.0;
		for (int i_level = 0; i_level <= searchLevel; i_level++)
		{
			pow_2_next_i_level *= 0.5;
			if (i_level == 0)
			{
				// Special case for level 0
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					CELL(i_level, i_dim) = 0.5;
					RATIO(i_level, i_dim) = 1.0;
					SLOPE(i_level, i_dim) = 0.0;
				}
			}
			else if (i_level == 1)
			{
				twoToPowerLevel *= 2;
				// Special case for level 1
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					site[i_dim] *= 2;
					if (site[i_dim] < 1)
					{
						CELL(i_level, i_dim) = 0.0;
						RATIO(i_level, i_dim) = 1 - site[i_dim];
						SLOPE(i_level, i_dim) = -twoToPowerLevel;
						if (i_level < searchLevel)
						{
							CELL(i_level + 1, i_dim) = 0.25;
						}
					}
					else
					{
						CELL(i_level, i_dim) = 1.0;
						site[i_dim] -= 1;
						RATIO(i_level, i_dim) = site[i_dim];
						SLOPE(i_level, i_dim) = twoToPowerLevel;
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = 0.75;
					}
				} // i_dim
			} // i_level==1
			else
			{
				twoToPowerLevel *= 2;
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					site[i_dim] *= 2;
					if (site[i_dim] < 1)
					{
						RATIO(i_level, i_dim) = site[i_dim];
						SLOPE(i_level, i_dim) = twoToPowerLevel;
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = CELL(i_level, i_dim) - pow_2_next_i_level;
					}
					else
					{
						site[i_dim] -= 1;
						RATIO(i_level, i_dim) = (1 - site[i_dim]);
						SLOPE(i_level, i_dim) = -twoToPowerLevel;
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = CELL(i_level, i_dim) + pow_2_next_i_level;
					}
				}
			} // i_level != 1
		}
	}


	inline
	void search_adept_no_slope(double* site, double* cell, double* ratio, double* slope, int searchLevel, int nDim)
	{
#define CELL(m,n) cell[(m)*nDim+(n)]
#define RATIO(m,n) ratio[(m)*nDim+(n)]
#define SLOPE(m,n) slope[(m)*nDim+(n)]

		double pow_2_next_i_level = 1.0;
		double twoToPowerLevel = 1.0;
		for (int i_level = 0; i_level <= searchLevel; i_level++)
		{
			pow_2_next_i_level *= 0.5;
			if (i_level == 0)
			{
				// Special case for level 0
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					CELL(i_level, i_dim) = 0.5;
					RATIO(i_level, i_dim) = 1.0;
				}
			}
			else if (i_level == 1)
			{
				twoToPowerLevel *= 2;
				// Special case for level 1
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					site[i_dim] *= 2;
					if (site[i_dim] < 1)
					{
						CELL(i_level, i_dim) = 0.0;
						RATIO(i_level, i_dim) = 1 - site[i_dim];
						if (i_level < searchLevel)
						{
							CELL(i_level + 1, i_dim) = 0.25;
						}
					}
					else
					{
						CELL(i_level, i_dim) = 1.0;
						site[i_dim] -= 1;
						RATIO(i_level, i_dim) = site[i_dim];
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = 0.75;
					}
				} // i_dim
			} // i_level==1
			else
			{
				twoToPowerLevel *= 2;
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					site[i_dim] *= 2;
					if (site[i_dim] < 1)
					{
						RATIO(i_level, i_dim) = site[i_dim];
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = CELL(i_level, i_dim) - pow_2_next_i_level;
					}
					else
					{
						site[i_dim] -= 1;
						RATIO(i_level, i_dim) = (1 - site[i_dim]);
						if (i_level < searchLevel)
							CELL(i_level + 1, i_dim) = CELL(i_level, i_dim) + pow_2_next_i_level;
					}
				}
			} // i_level != 1
		}
	}


	inline
	void loop_all_combinations_accumulate_surplus_adept(int i_vec, int nDim,
		int evalMaxLevel, GridInfoMap& info, std::vector<int>& setOfLevelCombinations,
		double* cell, double* ratio, double* slope, double* p_totalSurplus, double* p_gradient)
	{
		*p_totalSurplus = 0;
		for (int i_dim = 0; i_dim < nDim; i_dim++)
		{
			p_gradient[i_dim] = 0.0;
		}

		if (evalMaxLevel < 0)
			return;

		double gradientProd[ASG_MAX_DIM];
		inarray currentCell = {};

		int ptr = 0;
		while (ptr < setOfLevelCombinations.size())
		{
			int* levels = setOfLevelCombinations.data() + ptr;

			int totalLevel = sum_vec(levels, nDim);

			// Skip if exceeding max level
			if (totalLevel > evalMaxLevel)
				continue;

			// Construct current cell
			for (int i_dim = 0; i_dim < nDim; i_dim++)
			{
				currentCell[i_dim] = CELL(levels[i_dim], i_dim);
			}

			auto infoPiece = info.find(currentCell);

			if (infoPiece != info.end())
			{
				double ratioProd = 1.0;
				for (int i_dim = 0; i_dim < nDim ; i_dim++)
				{
					gradientProd[i_dim] = SLOPE(levels[i_dim], i_dim);
				}
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					ratioProd *= RATIO(levels[i_dim], i_dim);
					for (int j_dim = 0; j_dim < nDim ; j_dim++)
					{
						if (j_dim != i_dim)
							gradientProd[j_dim] *= RATIO(levels[i_dim], i_dim);
					}
				}
				double surplus = infoPiece->second.surplus[i_vec];
				p_totalSurplus[0] += surplus*ratioProd;
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					p_gradient[i_dim] += surplus * gradientProd[i_dim];
				}
			}

			ptr += nDim;
		}
	}

    inline
	void loop_all_combinations_accumulate_surplus_adept_no_grad(int i_vec, int nDim,
		int evalMaxLevel, GridInfoMap& info, std::vector<int>& setOfLevelCombinations,
		double* cell, double* ratio, double* slope, double* p_totalSurplus)
	{
		*p_totalSurplus = 0;
        
		if (evalMaxLevel < 0)
			return;

		inarray currentCell = {};

		int ptr = 0;
		while (ptr < setOfLevelCombinations.size())
		{
			int* levels = setOfLevelCombinations.data() + ptr;

			int totalLevel = sum_vec(levels, nDim);

			// Skip if exceeding max level
			if (totalLevel > evalMaxLevel)
				continue;

			// Construct current cell
			for (int i_dim = 0; i_dim < nDim; i_dim++)
			{
				currentCell[i_dim] = CELL(levels[i_dim], i_dim);
			}

			auto infoPiece = info.find(currentCell);

			if (infoPiece != info.end())
			{
				double ratioProd = 1.0;
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					ratioProd *= RATIO(levels[i_dim], i_dim);
				}
				double surplus = infoPiece->second.surplus[i_vec];
				p_totalSurplus[0] += surplus*ratioProd;
			}

			ptr += nDim;
		}
	}
    
	inline
	void loop_all_combinations_accumulate_surplus_vec_adept_no_grad(int nVec, int nDim,
		int evalMaxLevel, GridInfoMap& info, std::vector<int>& setOfLevelCombinations,
		double* cell, double* ratio, double* slope, double* p_totalSurplus)
	{
		for (int i_vec = 0; i_vec < nVec ; i_vec++)
		{
			p_totalSurplus[i_vec] = 0;
		}

		if (evalMaxLevel < 0)
			return;

		inarray currentCell = {};

		int ptr = 0;
		while (ptr < setOfLevelCombinations.size())
		{
			int* levels = setOfLevelCombinations.data() + ptr;

			int totalLevel = sum_vec(levels, nDim);

			// Skip if exceeding max level
			if (totalLevel > evalMaxLevel)
				continue;

			// Construct current cell
			for (int i_dim = 0; i_dim < nDim; i_dim++)
			{
				currentCell[i_dim] = CELL(levels[i_dim], i_dim);
			}

			auto infoPiece = info.find(currentCell);

			if (infoPiece != info.end())
			{
				double ratioProd = 1.0;
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					ratioProd *= RATIO(levels[i_dim], i_dim);
				}
				for (int i_vec = 0; i_vec < nVec; i_vec++)
				{
					double surplus = infoPiece->second.surplus[i_vec];
					p_totalSurplus[i_vec] += surplus*ratioProd;
				}
			}

			ptr += nDim;
		}
	}

	inline
	void loop_all_combinations_accumulate_surplus_vec_adept(int nVec, int nDim,
		int evalMaxLevel, GridInfoMap& info, std::vector<int>& setOfLevelCombinations,
		double* cell, double* ratio, double* slope, double* p_totalSurplus, double* p_gradient)
	{
		for (int i_vec = 0; i_vec < nVec ; i_vec++)
		{
			p_totalSurplus[i_vec] = 0;
			for (int i_dim = 0; i_dim < nDim ; i_dim++)
			{
				p_gradient[i_vec*nDim + i_dim] = 0.0;
			}
		}

		if (evalMaxLevel < 0)
			return;

		double gradientProd[ASG_MAX_DIM];
		inarray currentCell = {};

		int ptr = 0;
		while (ptr < setOfLevelCombinations.size())
		{
			int* levels = setOfLevelCombinations.data() + ptr;

			int totalLevel = sum_vec(levels, nDim);

			// Skip if exceeding max level
			if (totalLevel > evalMaxLevel)
				continue;

			// Construct current cell
			for (int i_dim = 0; i_dim < nDim; i_dim++)
			{
				currentCell[i_dim] = CELL(levels[i_dim], i_dim);
			}

			auto infoPiece = info.find(currentCell);

			if (infoPiece != info.end())
			{
				double ratioProd = 1.0;
				for (int i_dim = 0; i_dim < nDim ; i_dim++)
				{
					gradientProd[i_dim] = SLOPE(levels[i_dim], i_dim);
				}
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					ratioProd *= RATIO(levels[i_dim], i_dim);
					// Accumulate slope for other dimension
					for (int j_dim = 0; j_dim < nDim ; j_dim++)
					{
						if (j_dim != i_dim)
							gradientProd[j_dim] *= RATIO(levels[i_dim], i_dim);
					}
				}
				for (int i_vec = 0; i_vec < nVec; i_vec++)
				{
					double surplus = infoPiece->second.surplus[i_vec];
					p_totalSurplus[i_vec] += surplus*ratioProd;
					for (int i_dim = 0; i_dim < nDim ; i_dim++)
					{
						p_gradient[i_vec*nDim + i_dim] += surplus * gradientProd[i_dim];
					}
				}
			}

			ptr += nDim;
		}
	}
    
    inline
	void eval(int i_vec, double* cell, double* ratio, double* slope, int evalLevel, int nDim, GridInfoMap& info, std::vector<int>& setOfLevelCombinations, 
		double* fval)
	{
		loop_all_combinations_accumulate_surplus_adept_no_grad(i_vec, nDim, evalLevel, info, setOfLevelCombinations, cell, ratio, slope, fval);
	}

	inline
	void eval_adept(int i_vec, double* cell, double* ratio, double* slope, int evalLevel, int nDim, GridInfoMap& info, std::vector<int>& setOfLevelCombinations, 
		double* fval, double* gradient)
	{
		loop_all_combinations_accumulate_surplus_adept(i_vec, nDim, evalLevel, info, setOfLevelCombinations, cell, ratio, slope, fval, gradient);
	}

	inline
	void eval_vec_adept(int nVec, double* cell, double* ratio, double* slope, int evalLevel, int nDim, GridInfoMap& info, std::vector<int>& setOfLevelCombinations, 
		double* fval, double* gradient)
	{
		loop_all_combinations_accumulate_surplus_vec_adept(nVec, nDim, evalLevel, info, setOfLevelCombinations, cell, ratio, slope, fval, gradient);
	}

	inline
	void eval_vec_adept_no_grad(int nVec, double* cell, double* ratio, double* slope, int evalLevel, int nDim, GridInfoMap& info, std::vector<int>& setOfLevelCombinations, 
		double* fval)
	{
		loop_all_combinations_accumulate_surplus_vec_adept_no_grad(nVec, nDim, evalLevel, info, setOfLevelCombinations, cell, ratio, slope, fval);
	}

	class AsgInterpArrayAdoubleEvaluator
	{
	public:
		AsgInterpArray* interpArray;
		int numDim;
		int numVec;
		int numArray;
		int currentLevel; // Assume all interps share the same current level
		double* stateMin;
		double* stateRange;
		AsgInterp* interps;
		AsgInterpArrayAdoubleEvaluator(AsgInterpArray* _interpArray)
		{
			interpArray = _interpArray;
			interps = interpArray->interps.data();
			numDim = interpArray->numDim;
			numVec = interpArray->numVec;
			numArray = interpArray->numArray;
			currentLevel = interps[0].currentLevel;
			stateMin = interpArray->stateMin;
			stateRange = interpArray->stateRange;
		}

		void eval_adouble(int i_array, int i_vec, adouble* site, double* cell, double* ratio, double* slope,
			adouble* rslt, double* gradient)
		{
			double siteScaled[ASG_MAX_DIM];
			double fval[ASG_MAX_NVEC];
			for (int i_dim = 0; i_dim < numDim ; i_dim++)
			{
				siteScaled[i_dim] = (value(site[i_dim]) - stateMin[i_dim]) / stateRange[i_dim];
			}
			// search and eval
			search_adept_with_slope(siteScaled, cell, ratio, slope, currentLevel, numDim);
			AsgInterp& interp = interps[i_array];
			eval_adept(i_vec, cell, ratio, slope, currentLevel, numDim, interp.info, interp.setOfLevelCombinations, fval, gradient);

			rslt[0].set_value(fval[0]);
			for (int i_dim = 0; i_dim < numDim ; i_dim++)
			{
				gradient[i_dim] /= stateRange[i_dim];
			}
			rslt[0].add_derivative_dependence(site, gradient, numDim);
		}
        
        void eval_double(int i_array, int i_vec, double* site, double* cell, double* ratio, double* slope,
			double* rslt)
		{
			double siteScaled[ASG_MAX_DIM];
			for (int i_dim = 0; i_dim < numDim ; i_dim++)
			{
				siteScaled[i_dim] = (value(site[i_dim]) - stateMin[i_dim]) / stateRange[i_dim];
			}
			// search and eval
			search_adept_no_slope(siteScaled, cell, ratio, slope, currentLevel, numDim);
			AsgInterp& interp = interps[i_array];
			eval(i_vec, cell, ratio, slope, currentLevel, numDim, interp.info, interp.setOfLevelCombinations, rslt);
		}


		void scale_and_search_adouble(adouble* site, adouble* siteScaledTo01, double* cell, adouble* ratio)
		{
			for (int i_dim = 0; i_dim < numDim ; i_dim++)
			{
				siteScaledTo01[i_dim] = (site[i_dim] - stateMin[i_dim]) / stateRange[i_dim];
			}
			// search and eval
			search_adept(siteScaledTo01, cell, ratio, currentLevel, numDim);
		}

		void eval_vec_no_search_adouble(int i_array, adouble* siteScaledTo01, 
			double* cell, double* ratio, double* slope, adouble* rslt, double* fval, double* gradient)
		{
			AsgInterp& interp = interps[i_array];
			eval_vec_adept(numVec, cell, ratio, slope, currentLevel, numDim, interp.info, interp.setOfLevelCombinations, fval, gradient);
		}
        
        void eval_vec(int i_array, double* site, double* cell, double* ratio, double* slope, 
			double* rslt)
		{
			double siteScaled[ASG_MAX_DIM];
			for (int i_dim = 0; i_dim < numDim ; i_dim++)
			{
				siteScaled[i_dim] = (site[i_dim] - stateMin[i_dim]) / stateRange[i_dim];
			}
			// search and eval
			AsgInterp& interp = interps[i_array];
			search_adept_no_slope(siteScaled, cell, ratio, slope, currentLevel, numDim);
			eval_vec_adept_no_grad(numVec, cell, ratio, slope, currentLevel, numDim, interp.info, interp.setOfLevelCombinations, rslt);
		}

		void eval_vec_adouble(int i_array, adouble* site, double* cell, double* ratio, double* slope, 
			adouble* rslt, double* gradient, int evalGradFlag)
		{
			double siteScaled[ASG_MAX_DIM];
			double fval[ASG_MAX_NVEC];
			for (int i_dim = 0; i_dim < numDim ; i_dim++)
			{
				siteScaled[i_dim] = (value(site[i_dim]) - stateMin[i_dim]) / stateRange[i_dim];
			}
			// search and eval
			AsgInterp& interp = interps[i_array];
			if (evalGradFlag == 1)
			{
				search_adept_with_slope(siteScaled, cell, ratio, slope, currentLevel, numDim);
				eval_vec_adept(numVec, cell, ratio, slope, currentLevel, numDim, interp.info, interp.setOfLevelCombinations, fval, gradient);
			}
			else
			{
				search_adept_no_slope(siteScaled, cell, ratio, slope, currentLevel, numDim);
				eval_vec_adept_no_grad(numVec, cell, ratio, slope, currentLevel, numDim, interp.info, interp.setOfLevelCombinations, fval);
			}

			for (int i_vec = 0; i_vec < numVec ; i_vec++)
			{
				rslt[i_vec].set_value(fval[i_vec]);
				for (int i_dim = 0; i_dim < numDim ; i_dim++)
				{
					gradient[i_vec*numDim + i_dim] /= stateRange[i_dim];
				}
				rslt[i_vec].add_derivative_dependence(site, gradient + i_vec*numDim, numDim);
			}
		}
	};
}
