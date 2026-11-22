#include "functions.h"

#include <fenv.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

static bool expand_residuals(double *family, size_t width, size_t count,
                             double *column)
{
    const double *sum_errors = family + (count - 1) * width;
    for (size_t j = 0; j < width; ++j) {
        column[0] = family[j];
        column[1] = sum_errors[j];
        bool nonzero = column[0] != 0.0 || column[1] != 0.0;
        for (size_t row = 2; row < count; ++row) {
            column[row] = family[(row - 1) * width + j];
            nonzero |= column[row] != 0.0;
        }
        if (nonzero) {
            if (!accsumk(column, count, family + j, count, width)) {
                return false;
            }
        } else {
            for (size_t row = 0; row < count; ++row) {
                family[row * width + j] = 0.0;
            }
        }
    }
    return true;
}

static double residual_bound(const double *p, size_t ncoeff, double x)
{
    double bound = fabs(p[0]);
    for (size_t j = 1; j < ncoeff; ++j) {
        bound = bound * x + fabs(p[j]);
    }
    return bound;
}

static bool stop_test(double value, double theta, double bound, double gamma)
{
    const double u = 0x1p-53;
    const double b = (1 + 2*u)*bound*(1 + gamma/(2*(1 - 2*u)))/(1 - 4*u);
    if (b == 0.0) {
        return true;
    }
    const double lhs = (theta + b)/(1 - 2*u);
    const double rhs = u * ufp_bits(value);
    return isfinite(lhs) && lhs < rhs;
}

iter_comphorner_result iter_comphorner_accsum(const double *p, size_t ncoeff,
                                            double x, size_t max_iter)
{
    if (!p || ncoeff == 0 || max_iter == 0 || !isfinite(x)) {
        return (iter_comphorner_result){ITER_COMPHORNER_INVALID_INPUT, NAN, 0};
    }
    bool all_zero = true;
    for (size_t j = 0; j < ncoeff; ++j) {
        if (!isfinite(p[j])) {
            return (iter_comphorner_result){ITER_COMPHORNER_INVALID_INPUT, NAN, 0};
        }
        all_zero &= p[j] == 0.0;
    }
    if (ncoeff == 1 || x == 0.0 || all_zero) {
        return (iter_comphorner_result){ITER_COMPHORNER_OK, p[ncoeff - 1], 1};
    }
    /* 2^(2*M)*u <= 1, where M = ceil(log2(max_iter + 2)). */
    if (max_iter > (UINT32_C(1) << 26) - 2) {
        return (iter_comphorner_result){
            ITER_COMPHORNER_ITERATION_LIMIT_TOO_LARGE, NAN, 0};
    }
    if (fegetround() != FE_TONEAREST) {
        return (iter_comphorner_result){
            ITER_COMPHORNER_UNSUPPORTED_ROUNDING, NAN, 0};
    }

    const size_t width = ncoeff - 1;
    const double u = 0x1p-53;
    const double scaled_degree = 4*u*(double)(ncoeff - 2);
    if (scaled_degree >= 1.0) {
        return (iter_comphorner_result){ITER_COMPHORNER_INVALID_INPUT, NAN, 0};
    }
    const double gamma = scaled_degree / (1 - scaled_degree);

    double *family = NULL;
    double *sums = NULL;
    size_t capacity = 0;
    iter_comphorner_result result = {ITER_COMPHORNER_MAX_ITERATIONS, NAN, 0};

    for (size_t k = 1; k <= max_iter; ++k) {
        const size_t count = k + 1;
        if (count > capacity) {
            capacity = capacity == 0 ? 2 : 2 * capacity;
            if (capacity > max_iter + 1) {
                capacity = max_iter + 1;
            }
            if (width > SIZE_MAX / sizeof(*family) / capacity ||
                capacity > SIZE_MAX / sizeof(*sums) / 2) {
                result.status = ITER_COMPHORNER_ALLOCATION_FAILED;
                break;
            }
            double *new_family = realloc(family, capacity * width * sizeof(*family));
            if (!new_family) {
                result.status = ITER_COMPHORNER_ALLOCATION_FAILED;
                break;
            }
            family = new_family;
            double *new_sums = realloc(sums, 2 * capacity * sizeof(*sums));
            if (!new_sums) {
                result.status = ITER_COMPHORNER_ALLOCATION_FAILED;
                break;
            }
            sums = new_sums;
        }

        double *column = sums + capacity;
        double *sum_errors = family + k * width;
        double h;
        if (k == 1) {
            h = eft_horner(p, ncoeff, x, family, sum_errors);
        } else {
            /* Overwrite the first residual row as it is consumed. The new
             * residuals have one fewer coefficient and are right aligned. */
            h = eft_horner(family, width, x, family + 1, sum_errors + 1);
            family[0] = 0.0;
            sum_errors[0] = 0.0;
        }
        result.iterations = k;
        sums[k - 1] = h;
        if (!isfinite(h) || !expand_residuals(family, width, count, column)) {
            result.status = ITER_COMPHORNER_NONFINITE_INTERMEDIATE;
            break;
        }

        double value = h;
        double theta = 0.0;
        if (k == 2) {
            double error;
            two_sum(sums[0], sums[1], &value, &error);
            theta = fabs(error);
        } else if (k > 2) {
            /* Only the first two AccSumK terms enter the stopping bound. */
            double expansion[2];
            memcpy(column, sums, k * sizeof(*column));
            if (!accsumk(column, k, expansion, 2, 1)) {
                result.status = ITER_COMPHORNER_NONFINITE_INTERMEDIATE;
                break;
            }
            value = expansion[0];
            theta = fabs(expansion[1]) * (1 + 2*u) / (1 - 2*u);
        }
        const double bound = residual_bound(family, width, fabs(x));
        if (!isfinite(value) || !isfinite(theta) || !isfinite(bound)) {
            result.status = ITER_COMPHORNER_NONFINITE_INTERMEDIATE;
            break;
        }
        if (stop_test(value, theta, bound, gamma)) {
            result.status = ITER_COMPHORNER_OK;
            result.value = value;
            break;
        }
    }

    free(sums);
    free(family);
    return result;
}
