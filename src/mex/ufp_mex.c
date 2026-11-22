#include "mex.h"
#include "functions.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const double *x;
    double *y;
    mwSize count, i;
    mwSize ndim;
    mwSize *dims;
    mxArray *output;

    if (nrhs != 1) {
        mexErrMsgIdAndTxt("ufp_mex:InputCount", "Exactly one input is required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("ufp_mex:OutputCount", "At most one output is supported.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgIdAndTxt("ufp_mex:InputType", "Input must be a full real double array.");
    }
    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims = (mwSize *)mxMalloc(ndim * sizeof(*dims));
    memcpy(dims, mxGetDimensions(prhs[0]), ndim * sizeof(*dims));
    output = mxCreateUninitNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    mxFree(dims);
    x = mxGetDoubles(prhs[0]);
    y = mxGetDoubles(output);
    count = mxGetNumberOfElements(prhs[0]);
    for (i = 0; i < count; ++i) {
        y[i] = ufp_bits(x[i]);
    }
    if (nlhs == 0) {
        mxDestroyArray(output);
    } else {
        plhs[0] = output;
    }
}
