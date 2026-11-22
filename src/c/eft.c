#include "functions.h"

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

void two_sum(double a, double b, double *sum, double *error)
{
    const double s = a + b;
    const double z = s - a;
    *sum = s;
    *error = (a - (s - z)) + (b - z);
}

static void split_fp(double a, double *hi, double *lo)
{
    const double scaled = (0x1p27 + 1.0) * a;
    *hi = scaled - (scaled - a);
    *lo = a - *hi;
}

static void two_product(double a, double b, double *product, double *error)
{
    double ah, al, bh, bl;
    *product = a * b;
    split_fp(a, &ah, &al);
    split_fp(b, &bh, &bl);
    *error = al * bl - (((*product - ah * bh) - al * bh) - ah * bl);
}

double eft_horner(const double *p, size_t ncoeff, double x,
                  double *product_errors, double *sum_errors)
{
    double h = p[0];
    for (size_t j = 1; j < ncoeff; ++j) {
        double product, product_error, sum_error;
        two_product(h, x, &product, &product_error);
        two_sum(product, p[j], &h, &sum_error);
        product_errors[j - 1] = product_error;
        sum_errors[j - 1] = sum_error;
    }
    return h;
}
