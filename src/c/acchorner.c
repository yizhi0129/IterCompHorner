#define ITERCOMPHORNER_WITH_MPFR 1
#include "functions.h"

void acchorner_eval(mpfr_t y, const double *p, size_t ncoeff, double x)
{
    mpfr_t xx;
    mpfr_init2(xx, mpfr_get_prec(y));
    mpfr_set_d(xx, x, MPFR_RNDN);

    mpfr_set_d(y, p[0], MPFR_RNDN);
    for (size_t k = 1; k < ncoeff; ++k) {
        mpfr_mul(y, y, xx, MPFR_RNDN);
        mpfr_add_d(y, y, p[k], MPFR_RNDN);
    }

    mpfr_clear(xx);
}
