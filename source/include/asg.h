/**
	Adaptive Sparse Grid Interpolation with hat base function
	@author: Wenlan Luo, luowenlan@gmail.com
*/

#pragma once

#include <stdexcept>
#include "flat_hash_map.hpp"
#include <array>
#include <algorithm>
#include <vector>
#include <set>
#include <cmath>

#ifdef USE_OMP
#include "omp.h"
#endif

#ifndef ASG_MAX_DIM
#define ASG_MAX_DIM 10
#endif 

#ifndef ASG_MAX_NVEC
#define ASG_MAX_NVEC 100
#endif

#ifndef ASG_MAX_LEVEL
#define ASG_MAX_LEVEL 16
#endif

#define MAX(a,b)(((a)>(b))?(a):(b))

namespace AdaptiveSparseGridInterp {
	typedef std::array<int, ASG_MAX_DIM> levelarray;
	typedef std::array<double, ASG_MAX_DIM> inarray;
	typedef std::array<double, ASG_MAX_NVEC> outarray;

	struct GridInfo {
		levelarray levels = {};
		outarray surplus = {};
	};

	template <typename T>
	class CompareFirstDigits {
	private:
		int digits;
	public:
		CompareFirstDigits(int _digits) {
			digits = _digits;
		}

		bool operator()(const std::array<T, ASG_MAX_DIM>& a, const std::array<T, ASG_MAX_DIM>& b) const {
			for (int i = 0; i < digits ; i++)
			{
				if (a[i] < b[i])
					return true;
				else if (a[i] > b[i])
					return false;
			}
			return false;
		}
	};

	template <typename T>
	class EqualFirstDigits {
	private:
		int digits;
	public:
		EqualFirstDigits(int _digits)
		{
			digits = _digits;
		}

		constexpr bool operator()(const std::array<T, ASG_MAX_DIM>& a, const std::array<T, ASG_MAX_DIM>& b) const
		{
			for (int i = 0; i < digits ; i++)
			{
				if (a[i] != b[i])
					return false;
			}
			return true;
		}
	};

	//pow
	template<std::size_t n>
	struct helper_pow {
		inline static int pow2() {
			return 2 * helper_pow<n - 1>::pow2();
		}
	};

	//final specialization pow 
	template<>
	struct helper_pow<0> {
		inline static int pow2() {
			return 1;
		}
	};

	// a simple hash function that guarantees no collision in this special case
	template <typename T>
	class HashFirstDigits {
	private:
		int digits;
	public:
		HashFirstDigits(int _digits) {
			digits = _digits;
		}
		uint64_t operator()(const std::array<T, ASG_MAX_DIM>& a) const
		{
			uint64_t hash = 1469598103934665603ULL; // FNV offset
			for (int d = 0; d < digits; d++) {
				// Multiply coordinate[d] by 65536, cast to 32-bit
				uint32_t fi = (uint32_t)(a[d] * (1ULL << (ASG_MAX_LEVEL + 1)));
				// FNV combine
				hash ^= (uint64_t)fi;
				hash *= 1099511628211ULL;  // FNV prime
			}
			return hash;
		}
	};

	typedef ska::flat_hash_map<inarray, GridInfo, HashFirstDigits<double>, EqualFirstDigits<double>> GridInfoMap;
	typedef std::set<levelarray, CompareFirstDigits<int>> LevelSet;

	inline int find_pow_base_2(double gridPoint)
	{
		int powBase2 = 0;
		double diffGridPointToInt = 1;
		do
		{
			gridPoint *= 2;;
			diffGridPointToInt = gridPoint - (int)gridPoint;
			powBase2++;
		} while (diffGridPointToInt > 0);

		return powBase2;
	}

	inline int calculate_level(double gridPoint)
	{
		if (gridPoint == 0.5)
			return 0;
		else if (gridPoint == 0.0 || gridPoint == 1.0)
			return 1;
		else
		{
			return find_pow_base_2(gridPoint);
		}
	}

	inline int sum_vec(const int* vec, int dim)
	{
		int s = 0;
		for (int i = 0; i < dim; i++)
		{
			s += vec[i];
		}
		return s;
	}

	inline double max_abs(double* vec, int dim)
	{
		double s = 0;
		for (int i = 0; i < dim; i++)
		{
			s = MAX(fabs(vec[i]), s);
		}
		return s;
	}

	void loop_all_combinations_accumulate_surplus_vec(int nVec, int nDim,
		int evalMaxLevel, GridInfoMap& info, std::vector<int>& setOfLevelCombinations,
		inarray* cell, inarray* ratio, double* totalSurplus)
	{
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
				currentCell[i_dim] = cell[levels[i_dim]][i_dim];
			}

			auto infoPiece = info.find(currentCell);

			if (infoPiece != info.end())
			{
				double ratioProd = 1.0;
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					ratioProd *= ratio[levels[i_dim]][i_dim];
				}
				for (int i_vec = 0; i_vec < nVec; i_vec++)
				{
					totalSurplus[i_vec] += infoPiece->second.surplus[i_vec] * ratioProd;
				}
			}

			ptr += nDim;
		}

	}

	double loop_all_combinations_accumulate_surplus(int i_vec, int nDim,
		int evalMaxLevel, GridInfoMap& info, std::vector<int>& setOfLevelCombinations,
		inarray* cell, inarray* ratio)
	{
		double totalSurplus = 0;

		if (evalMaxLevel < 0)
			return totalSurplus;

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
				currentCell[i_dim] = cell[levels[i_dim]][i_dim];
			}

			auto infoPiece = info.find(currentCell);

			if (infoPiece != info.end())
			{
				double ratioProd = 1.0;
				for (int i_dim = 0; i_dim < nDim; i_dim++)
				{
					ratioProd *= ratio[levels[i_dim]][i_dim];
				}
				totalSurplus += infoPiece->second.surplus[i_vec] * ratioProd;
			}

			ptr += nDim;
		}

		return totalSurplus;
	}

	class AsgInterp {
	public:
		int numDim;
		int numVec;
		int currentLevel;

		GridInfoMap info = GridInfoMap(16, HashFirstDigits<double>(numDim), EqualFirstDigits<double>(numDim));
		std::vector<inarray> gridsCurrentLevel;
		std::vector<inarray> gridsNextLevel;

		std::vector<int> setOfLevelCombinations;

		void eval_vec(inarray* cell, inarray* ratio, int evalLevel, double* totalSurplus)
		{
			// evaluate recursively and accumulate surpluses
			memset(totalSurplus, 0, sizeof(double)*numVec);
			loop_all_combinations_accumulate_surplus_vec(numVec, numDim, evalLevel, info, setOfLevelCombinations, cell, ratio, totalSurplus);
		}

		double eval(int i_vec, inarray* cell, inarray* ratio, int evalLevel)
		{
			// evaluate recursively and accumulate surpluses
			return loop_all_combinations_accumulate_surplus(i_vec, numDim, evalLevel, info, setOfLevelCombinations, cell, ratio);
		}

		void search(inarray& site, inarray* cell, inarray* ratio, int searchLevel)
		{
			double pow_2_next_i_level = 1.0;
			for (int i_level = 0; i_level <= searchLevel; i_level++)
			{
				pow_2_next_i_level *= 0.5;
				if (i_level == 0)
				{
					// Special case for level 0
					for (int i_dim = 0; i_dim < numDim; i_dim++)
					{
						cell[i_level][i_dim] = 0.5;
						ratio[i_level][i_dim] = 1.0;
					}
				}
				else if (i_level == 1)
				{
					// Special case for level 1
					for (int i_dim = 0; i_dim < numDim; i_dim++)
					{
						site[i_dim] *= 2;
						if (site[i_dim] < 1)
						{
							cell[i_level][i_dim] = 0.0;
							ratio[i_level][i_dim] = 1 - site[i_dim];
							cell[i_level + 1][i_dim] = 0.25;
						}
						else
						{
							cell[i_level][i_dim] = 1.0;
							site[i_dim] -= 1;
							ratio[i_level][i_dim] = site[i_dim];
							cell[i_level + 1][i_dim] = 0.75;
						}
					} // i_dim
				} // i_level==1
				else
				{
					for (int i_dim = 0; i_dim < numDim; i_dim++)
					{
						site[i_dim] *= 2;
						if (site[i_dim] < 1)
						{
							ratio[i_level][i_dim] = site[i_dim];
							cell[i_level + 1][i_dim] = cell[i_level][i_dim] - pow_2_next_i_level;
						}
						else
						{
							site[i_dim] -= 1;
							ratio[i_level][i_dim] = (1 - site[i_dim]);
							cell[i_level + 1][i_dim] = cell[i_level][i_dim] + pow_2_next_i_level;
						}
					}
				} // i_level != 1
			}
		}

		inline double search_and_eval_at_level(int i_vec, double* _site, int evalLevel)
		{
			inarray site = {};
			memcpy(site.data(), _site, numDim * sizeof(double));

			// Prepare space and search
			inarray cell[ASG_MAX_LEVEL + 2];
			inarray ratio[ASG_MAX_LEVEL + 2];

			search(site, cell, ratio, evalLevel);
			return eval(i_vec, cell, ratio, evalLevel);
		}

		double search_and_eval(int i_vec, double* _site)
		{
			return search_and_eval_at_level(i_vec, _site, currentLevel);
		}

		inline void search_and_eval_vec_at_level(double* _site, double* rslt, int evalLevel)
		{
			inarray site = {};
			memcpy(site.data(), _site, numDim * sizeof(double));

			// Prepare space and search
			inarray cell[ASG_MAX_LEVEL + 2];
			inarray ratio[ASG_MAX_LEVEL + 2];

			search(site, cell, ratio, evalLevel);
			eval_vec(cell, ratio, evalLevel, rslt);
		}

		void search_and_eval_vec(double* _site, double* rslt)
		{
			search_and_eval_vec_at_level(_site, rslt, currentLevel);
		}

		void search_and_eval_batch_vec_omp(inarray* batchSites, int nSites, int evalLevel, double* rslt)
		{
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i_site = 0; i_site < nSites; i_site++)
			{
				inarray site = batchSites[i_site];
				inarray cell[ASG_MAX_LEVEL + 2];
				inarray ratio[ASG_MAX_LEVEL + 2];

				search(site, cell, ratio, evalLevel);
				eval_vec(cell, ratio, evalLevel, rslt + i_site*numVec);
			}
		}

		void search_and_eval_batch_vec_raw_array_omp(double* batchSites, int nSites, int evalLevel, double* rslt)
		{
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i_site = 0; i_site < nSites; i_site++)
			{
				inarray site;
				memcpy(site.data(), batchSites + i_site*numDim, numDim * sizeof(double));
				inarray cell[ASG_MAX_LEVEL + 2];
				inarray ratio[ASG_MAX_LEVEL + 2];

				search(site, cell, ratio, evalLevel);
				eval_vec(cell, ratio, evalLevel, rslt + i_site*numVec);
			}
		}

		void update_level_combinations()
		{
			// Call this when info changes
			LevelSet setOfLevelCombinationsSet = LevelSet(CompareFirstDigits<int>(numDim));
			setOfLevelCombinationsSet.clear();
			for (GridInfoMap::iterator it = info.begin(); it != info.end(); ++it)
				setOfLevelCombinationsSet.insert(it->second.levels);

			setOfLevelCombinations.clear();
			setOfLevelCombinations.reserve(numDim*setOfLevelCombinationsSet.size());
			for (auto it = setOfLevelCombinationsSet.begin(); it != setOfLevelCombinationsSet.end(); ++it)
			{
				setOfLevelCombinations.insert(setOfLevelCombinations.end(), it->begin(), it->begin() + numDim);
			}
		}

		void push_to_info(double* _fval)
		{
			std::vector<double> rslt;
			rslt.reserve(gridsNextLevel.size()*numVec);
			double* _rslt = rslt.data();
			search_and_eval_batch_vec_omp(gridsNextLevel.data(), gridsNextLevel.size(), currentLevel, _rslt);

			// Calculate surplus using vectorization
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < gridsNextLevel.size()*numVec; i++)
			{
				_rslt[i] = _fval[i] - _rslt[i];
			}

			// push to surpluses
			// auto it = surplus.begin();
			for (int i_grid = 0; i_grid < gridsNextLevel.size(); i_grid++)
			{
				// Calculate levels of each grid
				GridInfo infoAtCurrentGrid;
				const auto currentGrid = gridsNextLevel[i_grid];
				for (int i_dim = 0; i_dim < numDim ; i_dim++)
				{
					infoAtCurrentGrid.levels[i_dim] = calculate_level(currentGrid[i_dim]);
				}

				// Copy surplus
				memcpy(infoAtCurrentGrid.surplus.data(), _rslt + i_grid*numVec, numVec * sizeof(double));

#ifdef _WIN32
				info.emplace(gridsNextLevel[i_grid], infoAtCurrentGrid);
#else
				info[gridsNextLevel[i_grid]] = infoAtCurrentGrid;
#endif
			}

			currentLevel++;
			gridsCurrentLevel = gridsNextLevel;

			update_level_combinations();
		}

		void push_to_info_at_valid(double* _fval, bool* _valid)
		{
			std::vector<double> rslt;
			rslt.reserve(gridsNextLevel.size()*numVec);
			double* _rslt = rslt.data();
			search_and_eval_batch_vec_omp(gridsNextLevel.data(), gridsNextLevel.size(), currentLevel, _rslt);

			// Calculate surplus using vectorization
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < gridsNextLevel.size()*numVec; i++)
			{
				_rslt[i] = _fval[i] - _rslt[i];
			}

			// push to surpluses
			gridsCurrentLevel.clear();
			for (int i_grid = 0; i_grid < gridsNextLevel.size(); i_grid++)
			{
				if (!_valid[i_grid])
					continue;

				// Calculate levels of each grid
				GridInfo infoAtCurrentGrid;
				const auto currentGrid = gridsNextLevel[i_grid];
				for (int i_dim = 0; i_dim < numDim ; i_dim++)
				{
					infoAtCurrentGrid.levels[i_dim] = calculate_level(currentGrid[i_dim]);
				}

				// Copy surplus
				memcpy(infoAtCurrentGrid.surplus.data(), _rslt + i_grid*numVec, numVec * sizeof(double));

#ifdef _WIN32
				info.emplace(gridsNextLevel[i_grid], infoAtCurrentGrid);
#else
				info[gridsNextLevel[i_grid]] = infoAtCurrentGrid;
#endif
				gridsCurrentLevel.push_back(currentGrid);
			}

			currentLevel++;

			update_level_combinations();
		}

		void push_info_to_grids(double* _grids, int numGrids, double* _fval, int _currentLevel)
		{
			// Overwrite gridsNextLevel with grids and push _fval to info
			gridsNextLevel.clear();
			gridsNextLevel.reserve(numGrids);
			for (int i_grid = 0; i_grid < numGrids ; i_grid++)
			{
				inarray grid = {};
				memcpy(grid.data(), _grids + i_grid*numDim, sizeof(double)*numDim);
				gridsNextLevel.push_back(grid);
			}

			// reset currentLevel
			currentLevel = _currentLevel;
			push_to_info(_fval);
		}

		void construct_grids_next_level(double threshhold)
		{
			gridsNextLevel.clear();

			// Find the neighbor points and push to gridsNextLevel
			if (currentLevel == -1)
			{
				inarray gridL0 = {};
				for (int i_dim = 0; i_dim < numDim; i_dim++)
				{
					gridL0[i_dim] = 0.5;
				}
				gridsNextLevel.push_back(gridL0);
			}
			else
			{
				for (int i_grid = 0; i_grid < gridsCurrentLevel.size(); i_grid++)
				{
					inarray currentGrid = gridsCurrentLevel[i_grid];
					// Short-circuiting, if threshhold is 0.0, then add grids unconditionally
					if (threshhold == 0.0 || max_abs(info.at(currentGrid).surplus.data(), numVec) >= threshhold)
					{
						for (int i_dim = 0; i_dim < numDim; i_dim++)
						{
							inarray nextGrid = currentGrid;
							if (currentGrid[i_dim] == 0.5)
							{
								nextGrid[i_dim] = 0.0;
								gridsNextLevel.push_back(nextGrid);
								nextGrid[i_dim] = 1.0;
								gridsNextLevel.push_back(nextGrid);
							}
							else if (currentGrid[i_dim] == 0.0)
							{
								nextGrid[i_dim] = 0.25;
								gridsNextLevel.push_back(nextGrid);
							}
							else if (currentGrid[i_dim] == 1.0)
							{
								nextGrid[i_dim] = 0.75;
								gridsNextLevel.push_back(nextGrid);
							}
							else
							{
								double neighborNextLevel = pow(0.5, find_pow_base_2(currentGrid[i_dim]) + 1);
								nextGrid[i_dim] = currentGrid[i_dim] - neighborNextLevel;
								gridsNextLevel.push_back(nextGrid);
								nextGrid[i_dim] = currentGrid[i_dim] + neighborNextLevel;
								gridsNextLevel.push_back(nextGrid);
							}
						}
					}
				}
			}
			// Remove duplicates
			std::sort(gridsNextLevel.begin(), gridsNextLevel.end());
			gridsNextLevel.erase(unique(gridsNextLevel.begin(), gridsNextLevel.end()), gridsNextLevel.end());
		}

		/// Construct grid and surplus one level more
		template<typename Lambda>
		void construct_grids_and_surplus_using_f_next_level(Lambda f, double threshhold)
		{
			construct_grids_next_level(threshhold);

			// Evaluate function values at the new grids
			std::vector<double> fval;
			fval.reserve(gridsNextLevel.size()*numVec);
			double* _fval = fval.data();
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i_grid = 0; i_grid < gridsNextLevel.size(); i_grid++)
			{
				f(gridsNextLevel[i_grid].data(), _fval + i_grid*numVec);
			}

			// Construct surplus
			push_to_info(_fval);
		}

		// Construct grid and surplus from level -1 to finalLevel
		template<typename Lambda>
		void construct_grids_and_surplus_using_f_until_level(Lambda f, int _finalLevel, double threshhold)
		{
			while (currentLevel < _finalLevel)
			{
				construct_grids_and_surplus_using_f_next_level(f, threshhold);
			}
		}

		AsgInterp(int _numDim, int _numVec) :numDim(_numDim), numVec(_numVec), currentLevel(-1) {}
	};

	class AsgInterpArray
	{
	public:
		std::vector<AsgInterp> interps;
		int numDim;
		int numVec;
		int numArray;
		double stateMin[ASG_MAX_DIM] = {};
		double stateRange[ASG_MAX_DIM] = {};

		AsgInterpArray(int _numDim, int _numVec, int _numArray, double* _stateMin, double* _stateRange) : numDim(_numDim), numVec(_numVec), numArray(_numArray)
		{
			memcpy(stateMin, _stateMin, sizeof(double)*numDim);
			memcpy(stateRange, _stateRange, sizeof(double)*numDim);
			interps = std::vector<AsgInterp>(_numArray, AsgInterp(_numDim, _numVec));
		}

		double eval(int i_array, int i_vec, double* site)
		{
			double siteScaledTo01[ASG_MAX_DIM] = {};
			for (int i_dim = 0; i_dim < numDim; i_dim++)
			{
				siteScaledTo01[i_dim] = (site[i_dim] - stateMin[i_dim]) / stateRange[i_dim];
			}
			return interps[i_array].search_and_eval(i_vec, siteScaledTo01);
		}

		void eval_vec(int i_array, double* site, double* rslt)
		{
			double siteScaledTo01[ASG_MAX_DIM] = {};
			for (int i_dim = 0; i_dim < numDim; i_dim++)
			{
				siteScaledTo01[i_dim] = (site[i_dim] - stateMin[i_dim]) / stateRange[i_dim];
			}
			interps[i_array].search_and_eval_vec(siteScaledTo01, rslt);
		}
	};
}

#undef MAX

