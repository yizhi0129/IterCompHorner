#include "functions.h"

#include <float.h>
#include <math.h>

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

static double extract_vector(double sigma, double *values, size_t n)
{
    double tau = 0.0;
    for (size_t j = 0; j < n; ++j) {
        const double q = (sigma + values[j]) - sigma;
        values[j] -= q;
        tau += q;
    }
    return tau;
}

/* Transform followed by the TransformK leading-term reconstruction. */
static bool transform_k(double *values, size_t n, int m, double phi,
                        double *rho, double *res)
{
    double mu = 0.0;
    for (size_t j = 0; j < n; ++j) {
        if (!isfinite(values[j])) {
            return false;
        }
        mu = fmax(mu, fabs(values[j]));
    }
    if (mu == 0.0) {
        *res = *rho + 0.0;
        *rho = 0.0;
        return true;
    }

    double sigma = ldexp(1.0, m + (int)ceil(log2(mu)));
    if (!isfinite(sigma)) {
        return false;
    }
    const double scale = ldexp(1.0, m - 53);
    double t = *rho;
    for (;;) {
        const double tau = extract_vector(sigma, values, n);
        const double next_t = t + tau;
        if (!isfinite(next_t)) {
            return false;
        }
        /* MATLAB's realmin * 2^-53 rounds to zero in binary64. */
        if (fabs(next_t) >= phi * sigma || sigma == 0.0) {
            double tau1, tau2, residual = 0.0;
            two_sum(t, tau, &tau1, &tau2);
            for (size_t j = 0; j < n; ++j) {
                residual += values[j];
            }
            *res = tau1 + (tau2 + residual);
            *rho = tau2 - (*res - tau1);
            return isfinite(*res) && isfinite(*rho);
        }
        t = next_t;
        sigma *= scale;
    }
}

bool accsumk(double *values, size_t n, double *expansion,
             size_t terms, size_t stride)
{
    const int m = (int)ceil(log2((double)n + 2.0));
    if (m > 52) {
        return false;
    }
    const double phi = ldexp(1.0, 2 * m - 53);
    double rho = 0.0;
    for (size_t k = 0; k < terms; ++k) {
        double *term = &expansion[k * stride];
        if (!transform_k(values, n, m, phi, &rho, term)) {
            return false;
        }
        if (fabs(*term) <= 2.0 * DBL_MIN) {
            for (size_t j = k + 1; j < terms; ++j) {
                expansion[j * stride] = 0.0;
            }
            break;
        }
    }
    return true;
}
