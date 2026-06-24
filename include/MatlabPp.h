/* 
Convert Matlab Pp to C Struct
@Author: Wenlan Luo
@Date: 9/29/2016
*/

#include "mex.h"
#include "interp_lite.h"

template <int xDim, int interpOrder>
class MatlabPp {
public:
	double* xGrid[xDim];
	int xPts[xDim];
	int yDim;
	double* coefs;
	int coefsSize[xDim + 1];

    template <class DBL>
	DBL eval_4(DBL* xSite, int* cellOfSite, int vecIdx)
	{
		return search_eval4(xDim, xPts, xGrid, coefs, coefsSize, xSite, cellOfSite, vecIdx);
	}

    template <class DBL>
	DBL eval_2(DBL* xSite, int* cellOfSite, int vecIdx)
	{
		return search_eval2(xDim, xPts, xGrid, coefs, coefsSize, xSite, cellOfSite, vecIdx);
	}

    template <class DBL>
	DBL nosearch_eval_2(DBL* xSiteToLeft, int* cellOfSite, int vecIdx)
	{
		return receval2(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

    template <class DBL>
	DBL nosearch_eval_4(DBL* xSiteToLeft, int* cellOfSite, int vecIdx)
	{
		return receval4(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

	void search(adouble* xSite, int* cellOfSite, adouble* xSiteToLeft)
	{
		for (int j = 0; j < xDim; ++j)
		{
			hunt(xGrid[j], xPts[j] - 2, xSite[j].value(), cellOfSite + j);
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
	}
    
    void search(double* xSite, int* cellOfSite, double* xSiteToLeft)
	{
		for (int j = 0; j < xDim; ++j)
		{
			hunt(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
	}

	MatlabPp(const mxArray* pp)
	{
		const mxArray* ppBreaks = mxGetField(pp, 0, "breaks");
		const mxArray* ppCoefs = mxGetField(pp, 0, "coefs");
		const mxArray* ppDim = mxGetField(pp, 0, "dim");

		coefs = mxGetPr(ppCoefs);

		yDim = *((int*)mxGetData(ppDim));

		for (int j = 0; j < xDim; j++)
		{
			const mxArray* mxCurrentBreaks = mxGetCell(ppBreaks, j);
			xGrid[j] = (double*)mxGetPr(mxCurrentBreaks);
			xPts[j] = mxGetNumberOfElements(mxCurrentBreaks);
		}

		coefsSize[0] = yDim;
		for (int i = 0; i < xDim; ++i) {
			coefsSize[i + 1] = (xPts[i] - 1) * interpOrder;
		}
	}
};
