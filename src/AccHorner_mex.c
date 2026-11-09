#include "mex.h"

#include <mpfr.h>
#include <stdbool.h>

static void require_real_double(const mxArray *value, const char *name)
{
    if (!mxIsDouble(value) || mxIsComplex(value)) {
        mexErrMsgIdAndTxt("AccHorner:InvalidInput",
                          "%s must be a real double array.", name);
    }
}

static mpfr_prec_t parse_precision(const mxArray *value)
{
    double precision_double = mxGetScalar(value);

    if (precision_double < MPFR_PREC_MIN || precision_double > MPFR_PREC_MAX ||
        precision_double != (double)((mpfr_prec_t)precision_double)) {
        mexErrMsgIdAndTxt("AccHorner:InvalidPrecision",
                          "precision must be an integer number of bits supported by MPFR.");
    }

    return (mpfr_prec_t)precision_double;
}

static void horner_scalar(mpfr_t y,
                          const double *p,
                          mwSize ncoeff,
                          double x,
                          mpfr_prec_t precision)
{
    mpfr_t xx;
    mpfr_init2(xx, precision);
    mpfr_set_d(xx, x, MPFR_RNDN);

    mpfr_set_d(y, p[0], MPFR_RNDN);
    for (mwSize k = 1; k < ncoeff; ++k) {
        mpfr_mul(y, y, xx, MPFR_RNDN);
        mpfr_add_d(y, y, p[k], MPFR_RNDN);
    }

    mpfr_clear(xx);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2 || nrhs > 3) {
        mexErrMsgIdAndTxt("AccHorner:WrongNumberOfInputs",
                          "AccHorner_mex requires p, x, and optional precision.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("AccHorner:WrongNumberOfOutputs",
                          "AccHorner_mex returns one output.");
    }

    require_real_double(prhs[0], "p");
    require_real_double(prhs[1], "x");

    if (mxGetNumberOfElements(prhs[0]) == 0) {
        mexErrMsgIdAndTxt("AccHorner:EmptyPolynomial",
                          "p must contain at least one coefficient.");
    }

    mpfr_prec_t precision = 106;
    if (nrhs == 3) {
        require_real_double(prhs[2], "precision");
        if (mxGetNumberOfElements(prhs[2]) != 1) {
            mexErrMsgIdAndTxt("AccHorner:InvalidPrecision",
                              "precision must be a scalar.");
        }
        precision = parse_precision(prhs[2]);
    }

    const double *p = mxGetDoubles(prhs[0]);
    const double *x = mxGetDoubles(prhs[1]);
    mwSize ncoeff = mxGetNumberOfElements(prhs[0]);
    mwSize nx = mxGetNumberOfElements(prhs[1]);

    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]),
                                   mxGetDimensions(prhs[1]),
                                   mxDOUBLE_CLASS,
                                   mxREAL);
    double *yout = mxGetDoubles(plhs[0]);

    mpfr_t y;
    mpfr_init2(y, precision);

    for (mwSize i = 0; i < nx; ++i) {
        horner_scalar(y, p, ncoeff, x[i], precision);
        yout[i] = mpfr_get_d(y, MPFR_RNDN);
    }

    mpfr_clear(y);
}
