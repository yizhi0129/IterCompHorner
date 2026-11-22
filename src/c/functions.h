#ifndef ITERCOMPHORNER_FUNCTIONS_H
#define ITERCOMPHORNER_FUNCTIONS_H

#include <stdbool.h>
#include <stddef.h>

/* Define ITERCOMPHORNER_WITH_MPFR before including this header to use the
 * fixed-precision MPFR Horner evaluator. */
#ifdef ITERCOMPHORNER_WITH_MPFR
#include <mpfr.h>

/* Evaluate a nonempty coefficient array in descending degree order.
 * The caller initializes y, whose precision sets the working precision. */
void acchorner_eval(mpfr_t y, const double *p, size_t ncoeff, double x);
#endif

/* Return the highest power of two <= abs(x) for IEEE 754 binary64 x.
 * Return positive zero for zero, positive infinity for infinity, and a
 * NaN with the original payload and the sign bit cleared for NaN. */
double ufp_bits(double x);

void two_sum(double a, double b, double *sum, double *error);

/* p is nonempty and ordered by descending degree. Each error array has
 * ncoeff - 1 entries and the two outputs must not overlap. An output that
 * overlaps p may start at p[0] or p[1]. */
double eft_horner(const double *p, size_t ncoeff, double x,
                  double *product_errors, double *sum_errors);

/* Sum n entries into the first terms of the AccSumK expansion.
 * Values are scratch storage. Output entries are stride doubles apart.
 * Require n >= 1, 1 <= terms <= n, and distinct input and output buffers.
 * Returns false if a nonfinite intermediate prevents the transform. */
bool accsumk(double *values, size_t n, double *expansion,
             size_t terms, size_t stride);

typedef enum {
    ITER_COMPHORNER_OK = 0,
    ITER_COMPHORNER_INVALID_INPUT,
    ITER_COMPHORNER_ITERATION_LIMIT_TOO_LARGE,
    ITER_COMPHORNER_UNSUPPORTED_ROUNDING,
    ITER_COMPHORNER_ALLOCATION_FAILED,
    ITER_COMPHORNER_NONFINITE_INTERMEDIATE,
    ITER_COMPHORNER_MAX_ITERATIONS
} iter_comphorner_status;

typedef struct {
    iter_comphorner_status status;
    double value;      /* NaN on failure. */
    size_t iterations; /* Iterations attempted. Exact shortcuts count as one. */
} iter_comphorner_result;

/* Evaluate p(x) with revised Algorithm 12, using coefficients in descending
 * degree order. Requires binary64 round-to-nearest with gradual underflow.
 * Build with strict floating point (no fast-math or implicit FMA).
 * The stopping bound inherits the MATLAB EFT/AccSumK assumptions, including
 * absence of overflow and validity of the EFT residuals under underflow. */
iter_comphorner_result iter_comphorner_accsum(const double *p, size_t ncoeff,
                                            double x, size_t max_iter);

#endif
