#include "mex.h"
#include <float.h>
#include <stdint.h>
#include <string.h>

#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024
#error "ufp_mex requires IEEE 754 binary64 doubles."
#endif

/* Return the highest power of two <= abs(x), using only integer operations.
 * memcpy avoids undefined behaviour from pointer type punning. */
static double ufp_bits(double x)
{
    const uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
    const uint64_t fraction_mask = UINT64_C(0x000fffffffffffff);
    uint64_t bits, exponent, fraction, result;
    double y;

    memcpy(&bits, &x, sizeof(bits));
    exponent = bits & exponent_mask;
    fraction = bits & fraction_mask;
    if (exponent == exponent_mask) {
        /* Preserve infinity or the NaN payload, and remove the sign. */
        result = exponent | fraction;
    } else if (exponent != 0) {
        /* The stored biased exponent already encodes 2^floor(log2(abs(x))). */
        result = exponent;
    } else if (fraction != 0) {
        /* Subnormals have no implicit leading 1: retain the top fraction bit. */
        result = UINT64_C(1);
        while ((fraction >>= 1) != 0) {
            result <<= 1;
        }
    } else {
        result = 0; /* Both signs of zero map to positive zero. */
    }
    memcpy(&y, &result, sizeof(y));
    return y;
}

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
    if (sizeof(double) != sizeof(uint64_t)) {
        mexErrMsgIdAndTxt("ufp_mex:Platform", "IEEE 754 binary64 doubles are required.");
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
