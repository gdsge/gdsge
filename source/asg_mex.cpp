/**
	Interface to wrap the C++ AdaptiveSparseGrid class
	@author: Wenlan Luo, luowenlan@gmail.com
*/

#ifndef ASG_MAX_DIM
#define ASG_MAX_DIM 100
#endif 

#ifndef ASG_MAX_NVEC
#define ASG_MAX_NVEC 100
#endif

#ifndef ASG_MAX_LEVEL
#define ASG_MAX_LEVEL 20
#endif

#include "mex.h"
#include "class_handle.hpp"

// #define __linux__
#include "asg.h"

using namespace AdaptiveSparseGridInterp;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");

	// Get constans
	if (!strcmp("get_mex_constants", cmd)) {
		if (nrhs != 1)
            mexErrMsgTxt("new: 1 input expected.");
		if (nlhs != 3)
            mexErrMsgTxt("new: 3 input expected.");

		plhs[0] = mxCreateDoubleScalar(ASG_MAX_DIM);
		plhs[1] = mxCreateDoubleScalar(ASG_MAX_NVEC);
		plhs[2] = mxCreateDoubleScalar(ASG_MAX_LEVEL);

		return;
	}

	// New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nrhs != 6)
            mexErrMsgTxt("new: 6 input expected.");
		// Get parameters
		int numDim = (int)(*mxGetPr(prhs[1]));
		int numVec = (int)(*mxGetPr(prhs[2]));
		int numArray = (int)(*mxGetPr(prhs[3]));
		double* _stateMin = mxGetPr(prhs[4]);
		double* _stateRange = mxGetPr(prhs[5]);

		if (numDim > ASG_MAX_DIM)
		{
			mexErrMsgTxt("numDim > ASG_MAX_DIM");
		}

		if (numVec > ASG_MAX_NVEC)
		{
			mexErrMsgTxt("numVec > ASG_MAX_NVEC");
		}

        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<AsgInterpArray>(new AsgInterpArray(numDim, numVec, numArray, _stateMin, _stateRange));
        return;
    }

	// Delete
	if (!strcmp("delete", cmd)) {
		// Destroy the C++ object
        destroyObject<AsgInterpArray>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
	}

	// Get the class instance pointer from the second input
    AsgInterpArray* interpArray = convertMat2Ptr<AsgInterpArray>(prhs[1]);
	int numArray = interpArray->numArray;
	int nDim = interpArray->numDim;
	int nVec = interpArray->numVec;
	double* stateMin = interpArray->stateMin;
	double* stateRange = interpArray->stateRange;

	// Some get method
	if (!strcmp("get_current_level", cmd)) {
		if (nlhs != 1 || nrhs != 2)
			mexErrMsgTxt("get_current_level: input and output arguments error");

		AsgInterp& interp = interpArray->interps[0];
		plhs[0] = mxCreateDoubleScalar((double)interp.currentLevel);
        return;
	}

	// Return the grids to be evaluated
	if (!strcmp("get_eval_grids", cmd)) {
        if (nrhs != 3) mexErrMsgTxt("get_eval_grids: 3 input expected.");
		double threshhold = *mxGetPr(prhs[2]);

		//
		AsgInterp& interp0 = interpArray->interps[0];
		if (interp0.currentLevel + 1 > ASG_MAX_LEVEL)
		{
			mexErrMsgTxt("proposed level > MAX_LEVEL");
		}

		// Initialize the counter
		plhs[1] = mxCreateDoubleMatrix(1, numArray, mxREAL);
		mxArray* evalGridsLength = plhs[1];
		double* _evalGridsLength = mxGetPr(evalGridsLength);

		std::vector<inarray> allGrids;
		for (int i_grid = 0; i_grid < numArray; i_grid++)
		{
			AsgInterp& interp = interpArray->interps[i_grid];
			interp.construct_grids_next_level(threshhold);
			std::vector<inarray>& grids = interp.gridsNextLevel;
			allGrids.insert(allGrids.end(), grids.begin(), grids.end());
			// Accumulate the counter
			_evalGridsLength[i_grid] = grids.size();
		}

		// Allocate space
		plhs[0] = mxCreateDoubleMatrix(nDim, allGrids.size(), mxREAL);
		mxArray* evalGrids = plhs[0];
		double* _evalGrids = mxGetPr(evalGrids);

		plhs[2] = mxCreateDoubleMatrix(nDim, allGrids.size(), mxREAL);
		mxArray* evalGridsUnscaled = plhs[2];
		double* _evalGridsUnscaled = mxGetPr(evalGridsUnscaled);
		for (int i_grid = 0; i_grid < allGrids.size() ; i_grid++)
		{
			// memcpy(_evalGrids + i*interpArray->numDim, allGrids[i].data(), interpArray->numDim * sizeof(double));
			double* unscaledGrids = allGrids[i_grid].data();
			for (int i_dim = 0; i_dim < nDim ; i_dim++)
			{
				_evalGrids[i_grid*nDim + i_dim] = unscaledGrids[i_dim] * stateRange[i_dim] + stateMin[i_dim];
				_evalGridsUnscaled[i_grid*nDim + i_dim] = unscaledGrids[i_dim];
			}
		}

		return;
	}

	// Push evaluation results
	if (!strcmp("push_eval_results", cmd))
	{
		if (nrhs != 3) mexErrMsgTxt("push_eval_results: 3 input expected.");
		if (nlhs != 0) mexErrMsgTxt("push_eval_results: 0 output expected.");

		// Get the input vector
		const mxArray* mx_results = prhs[2];
		int results_m = mxGetM(mx_results);
		int results_n = mxGetN(mx_results);
		double* _results = mxGetPr(mx_results);

		// Check if the size is correct
		// Count grids
		int totalNumGridPoints = 0;
		std::vector<int> numGridPoints;
		numGridPoints.reserve(numArray);
		for (int i = 0; i < numArray; i++)
		{
			AsgInterp& interp = interpArray->interps[i];
			int gridSize = interp.gridsNextLevel.size();
			totalNumGridPoints += gridSize;
			numGridPoints.push_back(gridSize);
		}

		// Check if num of rows is consistent
		if (results_m != nVec)
			mexErrMsgTxt("results do not agree with number of vector functions");
		if (results_n != totalNumGridPoints)
			mexErrMsgTxt("results do not agree with number of points to be evaluated");

		double* ptr = _results;
		for (int i = 0; i < numArray; i++)
		{
			AsgInterp& interp = interpArray->interps[i];
			interp.push_to_info(ptr);
			ptr += numGridPoints[i] * nVec;
		}
		return;
	}

	// Push evaluation results
	if (!strcmp("push_eval_results_at_valid", cmd))
	{
		if (nrhs != 4) mexErrMsgTxt("push_eval_results_at_valid: 4 input expected.");
		if (nlhs != 0) mexErrMsgTxt("push_eval_results_at_valid: 0 output expected.");

		// Get the input vector
		const mxArray* mx_results = prhs[2];
		int results_m = mxGetM(mx_results);
		int results_n = mxGetN(mx_results);
		double* _results = mxGetPr(mx_results);

		const mxArray* mx_valid = prhs[3];
		if (!mxIsLogical(mx_valid))
			mexErrMsgTxt("valid must be a logical array");
		int valid_m = mxGetM(mx_valid);
		int valid_n = mxGetN(mx_valid);
		bool* _valid = (bool*)mxGetData(mx_valid);

		// Check if the size is correct
		// Count grids
		int totalNumGridPoints = 0;
		std::vector<int> numGridPoints;
		numGridPoints.reserve(numArray);
		for (int i = 0; i < numArray; i++)
		{
			AsgInterp& interp = interpArray->interps[i];
			int gridSize = interp.gridsNextLevel.size();
			totalNumGridPoints += gridSize;
			numGridPoints.push_back(gridSize);
		}

		// Check if num of rows is consistent
		if (results_m != nVec)
			mexErrMsgTxt("results do not agree with number of vector functions");
		if (results_n != totalNumGridPoints)
			mexErrMsgTxt("results do not agree with number of points to be evaluated");
		if (valid_m != 1)
			mexErrMsgTxt("valid should be a row vector");
		if (valid_n != totalNumGridPoints)
			mexErrMsgTxt("valid do not agree with number of points to be evaluated");

		double* ptr = _results;
		bool* ptr_valid = _valid;
		for (int i = 0; i < numArray; i++)
		{
			AsgInterp& interp = interpArray->interps[i];
			interp.push_to_info_at_valid(ptr, ptr_valid);
			ptr += numGridPoints[i] * nVec;
			ptr_valid += numGridPoints[i];
		}
		return;
	}

	// Push evaluation results once for all
	if (!strcmp("push_eval_results_at_grids", cmd))
	{
		if (nrhs != 6) mexErrMsgTxt("push_eval_results_at_grids: 6 input expected.");
		if (nlhs != 0) mexErrMsgTxt("push_eval_results_at_grids: 0 output expected.");

		// Get the input vector
		const mxArray* mx_evalArrayIdx = prhs[2];
		int evalArrayIdx_m = mxGetM(mx_evalArrayIdx);
		int evalArrayIdx_n = mxGetN(mx_evalArrayIdx);
		double* _evalArrayIdx = mxGetPr(mx_evalArrayIdx);

		const mxArray* mx_evalGrids = prhs[3];
		int evalGrids_m = mxGetM(mx_evalGrids);
		int evalGrids_n = mxGetN(mx_evalGrids);
		double* _evalGrids = mxGetPr(mx_evalGrids);

		const mxArray* mx_results = prhs[4];
		int results_m = mxGetM(mx_results);
		int results_n = mxGetN(mx_results);
		double* _results = mxGetPr(mx_results);
		
		int currentLevel = (int)*mxGetPr(prhs[5]);

		// Check if the size is correct
		// Check if num of rows is consistent
		if (evalGrids_n != evalArrayIdx_n)
			mexErrMsgTxt("length of evalArrayIdx and evalGrids do not agree");
		if (evalArrayIdx_n != results_n)
			mexErrMsgTxt("length of evalArrayIdx and results do not agree");
		if (evalGrids_m != nDim)
			mexErrMsgTxt("evalGrids do not agree with nDim");
		if (results_m != nVec)
			mexErrMsgTxt("results do not agree with number of vector functions");

		// collect grids results for each array of interp
		std::vector<std::vector<double>> gridsAtArray(numArray);
		std::vector<std::vector<double>> resultsAtArray(numArray);
		for (int i_grid = 0; i_grid < evalGrids_n ; i_grid++)
		{
			int i_array = (int)(_evalArrayIdx[i_grid]) - 1;
			gridsAtArray[i_array].insert(gridsAtArray[i_array].end(), 
				_evalGrids + i_grid*nDim, _evalGrids + (i_grid + 1)*nDim);
			resultsAtArray[i_array].insert(resultsAtArray[i_array].end(), 
				_results + i_grid*nVec, _results + (i_grid + 1)*nVec);
		}
		
		// For each array, push grids and results
		for (int i_array = 0; i_array < numArray ; i_array++)
		{
			AsgInterp& interp = interpArray->interps[i_array];
			interp.push_info_to_grids(gridsAtArray[i_array].data(),
				gridsAtArray[i_array].size() / nDim, resultsAtArray[i_array].data(),
				currentLevel);
		}

		return;
	}

	// Return grids in a cell
	if (!strcmp("get_grids_info", cmd)) {
        if (nrhs != 2) mexErrMsgTxt("get_grids_info: 2 input expected.");
        if (nlhs != 3) mexErrMsgTxt("get_grids_info: 3 output expected.");

		// Construct the struct
		plhs[0] = mxCreateCellMatrix(1, numArray);
		plhs[1] = mxCreateCellMatrix(1, numArray);
		plhs[2] = mxCreateCellMatrix(1, numArray);
		mxArray* gridsCell = plhs[0];
		mxArray* surplusCell = plhs[1];
		mxArray* levelsCell = plhs[2];
		for (int i = 0; i < numArray ; i++)
		{
			// copy grids to mx_grid
			AsgInterp& interp = interpArray->interps[i];
			mxArray* mx_grid = mxCreateDoubleMatrix(nDim, interp.info.size(), mxREAL);
			mxArray* mx_surplus = mxCreateDoubleMatrix(nVec, interp.info.size(), mxREAL);
			mxArray* mx_levels = mxCreateDoubleMatrix(nDim, interp.info.size(), mxREAL);
			double* p_grid = mxGetPr(mx_grid);
			double* p_surplus = mxGetPr(mx_surplus);
			double* p_levels = mxGetPr(mx_levels);
			for (GridInfoMap::iterator it = interp.info.begin(); it != interp.info.end(); ++it)
			{
				memcpy(p_grid, it->first.data(), nDim * sizeof(double));
				auto info = it->second;
				memcpy(p_surplus, info.surplus.data(), nVec * sizeof(double));

				for (int i_dim = 0; i_dim < nDim ; i_dim++)
				{
					p_levels[i_dim] = (double)info.levels[i_dim];
				}

				p_grid += nDim;
				p_surplus += nVec;
				p_levels += nDim;
			}
			// set cell
			mxSetCell(gridsCell, i, mx_grid);
			mxSetCell(surplusCell, i, mx_surplus);
			mxSetCell(levelsCell, i, mx_levels);
		}
		return;
	}

	// Evaluate at site, index and array
	if (!strcmp("eval", cmd)) {
		if (nrhs != 5) mexErrMsgTxt("5 input expected.");
		if (nlhs != 1) mexErrMsgTxt("1 output expected.");

		// get input
		const mxArray* mx_arrayIdx = prhs[2];
		const mxArray* mx_vecIdx = prhs[3];
		const mxArray* mx_sites = prhs[4];

		// get data pointer
		double* _sites = mxGetPr(mx_sites);
		double* _vecIdx = mxGetPr(mx_vecIdx);
		double* _arrayIdx = mxGetPr(mx_arrayIdx);

		// get size of sites
		int sites_m = mxGetM(mx_sites);
		int numSites = mxGetN(mx_sites);

		if (numSites == 0 || sites_m == 0)
			mexErrMsgTxt("sites cannot be empty");
		if (sites_m != nDim)
			mexErrMsgTxt("sites dimension error");

		// Create output
		plhs[0] = mxCreateDoubleMatrix(1, numSites, mxREAL);
		double* _results = mxGetPr(plhs[0]);

#pragma omp parallel for
		for (int i = 0; i < numSites ; i++)
		{
			_results[i] = interpArray->eval((int)(_arrayIdx[i] - 1), (int)(_vecIdx[i] - 1), _sites + i*nDim);
		}

		return;
	}

	if (!strcmp("eval_vec", cmd)) {
		if (nrhs != 4) mexErrMsgTxt("4 input expected.");
		if (nlhs != 1) mexErrMsgTxt("1 output expected.");

		// get input
		const mxArray* mx_arrayIdx = prhs[2];
		const mxArray* mx_sites = prhs[3];

		// get data pointer
		double* _sites = mxGetPr(mx_sites);
		double* _arrayIdx = mxGetPr(mx_arrayIdx);

		// get size of sites
		int sites_m = mxGetM(mx_sites);
		int numSites = mxGetN(mx_sites);

		if (sites_m != nDim)
			mexErrMsgTxt("sites dimension error");

		// Create output
		plhs[0] = mxCreateDoubleMatrix(nVec, numSites, mxREAL);
		double* _results = mxGetPr(plhs[0]);

#pragma omp parallel for
		for (int i = 0; i < numSites ; i++)
		{
			interpArray->eval_vec((int)(_arrayIdx[i] - 1), _sites + i*nDim, _results + i*nVec);
		}
		return;
	}

	// Return current and next level grids information
	if (!strcmp("get_grids_current_and_next", cmd)) {
		if (nrhs != 2) mexErrMsgTxt("2 input expected.");
		if (nlhs != 2) mexErrMsgTxt("2 output expected.");

		plhs[0] = mxCreateCellMatrix(1, numArray);
		plhs[1] = mxCreateCellMatrix(1, numArray);

		for (int i_array = 0; i_array < numArray ; i_array++)
		{
			AsgInterp& interp = interpArray->interps[i_array];
			std::vector<inarray>& gridsCurrentLevel = interp.gridsCurrentLevel;
			std::vector<inarray>& gridsNextLevel = interp.gridsNextLevel;
			mxArray* mx_gridsCurrentLevel = mxCreateDoubleMatrix(nDim, gridsCurrentLevel.size(), mxREAL);
			double* _gridsCurrentLevel = mxGetPr(mx_gridsCurrentLevel);
			mxArray* mx_gridsNextLevel = mxCreateDoubleMatrix(nDim, gridsNextLevel.size(), mxREAL);
			double* _gridsNextLevel = mxGetPr(mx_gridsNextLevel);
			for (int i_grid = 0; i_grid < gridsCurrentLevel.size() ; i_grid++)
			{
				memcpy(_gridsCurrentLevel + i_grid*nDim, gridsCurrentLevel[i_grid].data(), sizeof(double)*nDim);
			}
			for (int i_grid = 0; i_grid < gridsNextLevel.size(); i_grid++)
			{
				memcpy(_gridsNextLevel + i_grid*nDim, gridsNextLevel[i_grid].data(), sizeof(double)*nDim);
			}
			// Set cell
			mxSetCell(plhs[0], i_array, mx_gridsCurrentLevel);
			mxSetCell(plhs[1], i_array, mx_gridsNextLevel);
		}
		return;
	}

	if (!strcmp("reconstruct_internal", cmd)) {
		if (nrhs != 8) mexErrMsgTxt("8 input expected.");
		if (nlhs != 0) mexErrMsgTxt("0 output expected.");

		// Get from input
		int currentLevel = (int)*mxGetPr(prhs[2]);
		const mxArray* mx_grids = prhs[3];
		const mxArray* mx_surplus = prhs[4];
		const mxArray* mx_levels = prhs[5];
		const mxArray* mx_gridsCurrentLevel = prhs[6];
		const mxArray* mx_gridsNextLevel = prhs[7];
		
		// The routine does not check consistency of input (i.e., grids should be cell), make sure this is only called from the matlab class file

		// Some size information
		int mx_grids_N = mxGetN(mx_grids);
		if (mx_grids_N != numArray)
			mexErrMsgTxt("input grids cell size does not agree with numArray");

		// Input information into interpArray
		for (int i_array = 0; i_array < numArray ; i_array++)
		{
			// Shared information
			AsgInterp& interp = interpArray->interps[i_array];
			interp.currentLevel = currentLevel;

			// Get from cell
			const mxArray* mx_grids_i_array = mxGetCell(mx_grids, i_array);
			int nGrids = mxGetN(mx_grids_i_array);
			double* _grids = mxGetPr(mx_grids_i_array);

			double* _surplus = mxGetPr(mxGetCell(mx_surplus, i_array));
			double* _levels = mxGetPr(mxGetCell(mx_levels, i_array));

			const mxArray* mx_gridsCurrentLevel_i_array = mxGetCell(mx_gridsCurrentLevel, i_array);
			int nGridsCurrentLevel = mxGetN(mx_gridsCurrentLevel_i_array);
			double* _gridsCurrentLevel = mxGetPr(mx_gridsCurrentLevel_i_array);
			
			const mxArray* mx_gridsNextLevel_i_array = mxGetCell(mx_gridsNextLevel, i_array);
			int nGridsNextLevel = mxGetN(mx_gridsNextLevel_i_array);
			double* _gridsNextLevel = mxGetPr(mx_gridsNextLevel_i_array);

			// Directly write to internal information. Unsafe but simple
			// grids
			// Get the begin
			auto it = interp.info.end();
			for (int i_grid = 0; i_grid < nGrids ; i_grid++)
			{
				inarray grid = {};
				memcpy(grid.data(), _grids + i_grid*nDim, sizeof(double)*nDim);
				GridInfo infoPiece;
				for (int i_dim = 0; i_dim < nDim ; i_dim++)
				{
					infoPiece.levels[i_dim] = (int)_levels[i_grid*nDim + i_dim];
				}
				memcpy(infoPiece.surplus.data(), _surplus + i_grid*nVec, sizeof(double)*nVec);

				it = interp.info.emplace_hint(it, grid, infoPiece);
			}
			interp.update_level_combinations();

			// Write the gridsCurrent and gridsNext information
			interp.gridsCurrentLevel.reserve(nGridsCurrentLevel);
			for (int i_grid = 0; i_grid < nGridsCurrentLevel ; i_grid++)
			{
				inarray grid = {};
				memcpy(grid.data(), _gridsCurrentLevel + i_grid*nDim, sizeof(double)*nDim);
				interp.gridsCurrentLevel.push_back(grid);
			}

			interp.gridsNextLevel.reserve(nGridsNextLevel);
			for (int i_grid = 0; i_grid < nGridsNextLevel ; i_grid++)
			{
				inarray grid = {};
				memcpy(grid.data(), _gridsNextLevel + i_grid*nDim, sizeof(double)*nDim);
				interp.gridsNextLevel.push_back(grid);
			}
		}
		return;
	}

	mexErrMsgTxt("Methods not implemented");
}

