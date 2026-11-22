#define _POSIX_C_SOURCE 200809L

#include <sollya.h>
#include "functions.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum {
    FIRST_EXPONENT = 2,
    LAST_EXPONENT = 35,
    MAX_DEGREE = 500,
    RESULT_PRECISION = 53,
    UNIT_ROUNDOFF_EXPONENT = 53,
    TRIALS = 10,
    MAX_ITERATIONS = 100,
    RANDOM_SEED = 2026
};

typedef struct {
    const double *p;
    size_t degree;
    double x;
    sollya_obj_t polynomial;
    mpfr_ptr point;
} evaluation_case;

typedef struct {
    double condition;
    double iter_seconds;
    double sollya_seconds;
    bool both_faithful;
    bool iter_error_lt_u;
    bool sollya_error_lt_u;
} trial_result;

static uint32_t next_random(uint32_t *state)
{
    uint32_t value = *state;
    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    *state = value;
    return value;
}

static void make_polynomial(double *p, size_t degree, double x, double target,
                            uint32_t *random_state)
{
    const double x_squared = x * x;

    /* Each block contributes (a*t + b)*(t^2 - x^2), times a power of t.
     * The point and coefficients have short binary significands, so the
     * products defining each block are exact. At t=x the blocks
     * cancel exactly and only the final linear term determines p(x). */
    for (size_t j = 0; j + 3 < degree - 1; j += 4) {
        const uint32_t draw = next_random(random_state);
        const double a = ldexp(1.0 + (double)(draw & 255U) / 256.0,
                               -(int)((draw >> 8) & 7U));
        const double b = ldexp(a, -50 - (int)((draw >> 11) & 7U));
        p[j] = draw & (1U << 14) ? a : -a;
        p[j + 1] = draw & (1U << 15) ? b : -b;
        p[j + 2] = -p[j] * x_squared;
        p[j + 3] = -p[j + 1] * x_squared;
    }

    double magnitude = 0.0;
    for (size_t j = 0; j <= degree; ++j) {
        magnitude = magnitude * x + fabs(p[j]);
    }
    /* cond(p,x) = (magnitude + tail*(x+1)) / (tail*(x+1)). */
    const double tail = magnitude / ((target - 1.0) * (x + 1.0));
    p[degree - 1] = tail;
    p[degree] = tail;
}

static void exact_reference(const double *p, size_t degree, double point,
                            mpq_t value, mpq_t magnitude)
{
    mpq_t x, coefficient;
    mpq_inits(x, coefficient, NULL);
    mpq_set_d(x, point);
    for (size_t j = 0; j <= degree; ++j) {
        mpq_set_d(coefficient, p[j]);
        mpq_mul(value, value, x);
        mpq_add(value, value, coefficient);
        mpq_abs(coefficient, coefficient);
        mpq_mul(magnitude, magnitude, x);
        mpq_add(magnitude, magnitude, coefficient);
    }
    mpq_clears(x, coefficient, NULL);
}

static sollya_obj_t make_sollya_polynomial(const double *p, size_t degree)
{
    sollya_obj_t polynomial = sollya_lib_constant_from_double(p[0]);
    /* The build_function calls take ownership of their operands. */
    for (size_t j = 1; j <= degree; ++j) {
        polynomial = sollya_lib_build_function_add(
            sollya_lib_build_function_mul(polynomial, sollya_lib_free_variable()),
            sollya_lib_constant_from_double(p[j]));
    }
    return polynomial;
}

static bool sollya_certified(sollya_fp_result_t status)
{
    return status & (SOLLYA_FP_FLAG_PROVEN_EXACT | SOLLYA_FP_FLAG_CORRECTLY_ROUNDED |
                     SOLLYA_FP_FLAG_FAITHFUL);
}

static double relative_error(mpfr_srcptr value, mpq_srcptr reference, mpq_t scratch,
                             bool *below_u)
{
    *below_u = false;
    if (!mpfr_number_p(value)) {
        return NAN;
    }
    mpfr_get_q(scratch, value);
    mpq_sub(scratch, scratch, reference);
    mpq_div(scratch, scratch, reference);
    mpq_abs(scratch, scratch);
    const double error = mpq_get_d(scratch);
    /* Compare error / u to one before any floating-point rounding. */
    mpq_mul_2exp(scratch, scratch, UNIT_ROUNDOFF_EXPONENT);
    *below_u = mpq_cmp_ui(scratch, 1, 1) < 0;
    return error;
}

static double measure_seconds(const evaluation_case *input, bool use_sollya,
                              mpfr_ptr sollya_result)
{
    /* Batch fast calls to keep individual clock reads out of the measurement. */
    for (size_t repetitions = 1; repetitions <= 1048576; repetitions *= 2) {
        struct timespec start, end;
        if (clock_gettime(CLOCK_MONOTONIC, &start) != 0) {
            return NAN;
        }
        if (use_sollya) {
            for (size_t i = 0; i < repetitions; ++i) {
                sollya_lib_evaluate_function_at_point(
                    sollya_result, input->polynomial, input->point, NULL);
            }
        } else {
            for (size_t i = 0; i < repetitions; ++i) {
                iter_comphorner_accsum(input->p, input->degree + 1, input->x,
                                       MAX_ITERATIONS);
            }
        }
        if (clock_gettime(CLOCK_MONOTONIC, &end) != 0) {
            return NAN;
        }
        const double elapsed = (double)(end.tv_sec - start.tv_sec) +
                               1e-9 * (double)(end.tv_nsec - start.tv_nsec);
        if (elapsed >= 0.001) {
            return elapsed / (double)repetitions;
        }
    }
    return NAN;
}

static bool run_trial(size_t degree, int exponent, size_t trial, uint32_t *random_state,
                      trial_result *sample, double *max_iter_error, double *max_sollya_error,
                      FILE *detail_output)
{
    const double target = pow(10.0, exponent);
    const double x = 1.0 + (double)(1U + (next_random(random_state) & 15U)) / 256.0;
    double p[MAX_DEGREE + 1] = {0.0};
    make_polynomial(p, degree, x, target, random_state);
    mpq_t reference, magnitude, ratio;
    mpq_inits(reference, magnitude, ratio, NULL);
    exact_reference(p, degree, x, reference, magnitude);
    if (mpq_sgn(reference) == 0) {
        fprintf(stderr, "Zero reference for degree %zu at 1e%d\n", degree, exponent);
        mpq_clears(reference, magnitude, ratio, NULL);
        return false;
    }
    mpq_abs(ratio, reference);
    mpq_div(ratio, magnitude, ratio);
    sample->condition = mpq_get_d(ratio);
    const bool condition_matches = fabs(sample->condition / target - 1.0) <= 1e-12;

    mpfr_t point, sollya_value, iter_value, downward, upward;
    mpfr_inits2(RESULT_PRECISION, point, sollya_value, iter_value, downward, upward,
                (mpfr_ptr)NULL);
    mpfr_set_d(point, x, MPFR_RNDN);
    sollya_obj_t polynomial = make_sollya_polynomial(p, degree);
    const evaluation_case input = {p, degree, x, polynomial, point};
    const sollya_fp_result_t sollya_status =
        sollya_lib_evaluate_function_at_point(sollya_value, polynomial, point, NULL);
    const iter_comphorner_result iter = iter_comphorner_accsum(p, degree + 1, x, MAX_ITERATIONS);
    mpfr_set_d(iter_value, iter.value, MPFR_RNDN);
    mpfr_set_q(downward, reference, MPFR_RNDD);
    mpfr_set_q(upward, reference, MPFR_RNDU);
    const bool iter_faithful = iter.status == ITER_COMPHORNER_OK &&
        (mpfr_equal_p(iter_value, downward) || mpfr_equal_p(iter_value, upward));
    const bool sollya_faithful = sollya_certified(sollya_status) &&
        (mpfr_equal_p(sollya_value, downward) || mpfr_equal_p(sollya_value, upward));
    sample->both_faithful = iter_faithful && sollya_faithful;
    const double iter_error = relative_error(iter_value, reference, ratio,
                                             &sample->iter_error_lt_u);
    const double sollya_error = relative_error(sollya_value, reference, ratio,
                                               &sample->sollya_error_lt_u);
    *max_iter_error = fmax(*max_iter_error, iter_error);
    *max_sollya_error = fmax(*max_sollya_error, sollya_error);

    /* Include time to return a failure. Alternate the method called first. */
    for (size_t order = 0; order < 2; ++order) {
        const bool use_sollya = (trial + order) % 2 != 0;
        const double seconds = measure_seconds(&input, use_sollya, sollya_value);
        if (use_sollya) {
            sample->sollya_seconds = seconds;
        } else {
            sample->iter_seconds = seconds;
        }
    }

    char row[512];
    const int row_length = snprintf(row, sizeof(row),
        "%zu\t1e%d\t%zu\t%.16e\t%.16e\t%.16e\t%d\t%zu\t%d\t%d\t%d\t%.9e\t%.9e\t%d\t%d\n",
        degree, exponent, trial + 1, sample->condition, iter_error,
        sollya_error, (int)iter.status, iter.iterations, (int)sollya_status,
        iter_faithful, sollya_faithful, sample->iter_seconds, sample->sollya_seconds,
        sample->iter_error_lt_u, sample->sollya_error_lt_u);
    bool row_written = row_length >= 0 && (size_t)row_length < sizeof(row);
    if (row_written) {
        fputs(row, stdout);
        row_written = fputs(row, detail_output) != EOF;
    }

    sollya_lib_clear_obj(polynomial);
    mpfr_clears(point, sollya_value, iter_value, downward, upward, (mpfr_ptr)NULL);
    mpq_clears(reference, magnitude, ratio, NULL);
    return row_written && condition_matches && sollya_faithful &&
           (iter_faithful || iter.status == ITER_COMPHORNER_MAX_ITERATIONS) &&
           isfinite(sample->iter_seconds) && isfinite(sample->sollya_seconds);
}

static int compare_doubles(const void *lhs, const void *rhs)
{
    const double a = *(const double *)lhs;
    const double b = *(const double *)rhs;
    return (a > b) - (a < b);
}

static double median(double samples[TRIALS])
{
    qsort(samples, TRIALS, sizeof(samples[0]), compare_doubles);
    return (samples[TRIALS / 2 - 1] + samples[TRIALS / 2]) / 2.0;
}

static bool close_report(FILE *report, const char *path)
{
    const bool write_failed = ferror(report) != 0;
    if (fclose(report) == EOF || write_failed) {
        fprintf(stderr, "Could not write %s\n", path);
        return false;
    }
    return true;
}

int main(int argc, char **argv)
{
    const size_t degrees[] = {5, 10, 50, 100, 300, MAX_DEGREE};
    const char *output_path = argc >= 2 ? argv[1] :
        "results/iter_comphorner_sollya_performance_times.txt";
    const char *detail_path = argc >= 4 ? argv[3] :
        "results/iter_comphorner_sollya_trial_details.txt";
    if (argc > 4) {
        fprintf(stderr, "Usage: %s [output.txt [seed [trial_details.txt]]]\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (strcmp(output_path, detail_path) == 0) {
        fputs("Timing and trial detail reports require distinct output paths.\n", stderr);
        return EXIT_FAILURE;
    }
    uint32_t seed = RANDOM_SEED;
    if (argc >= 3) {
        char *end;
        errno = 0;
        const uintmax_t parsed = strtoumax(argv[2], &end, 10);
        if (argv[2][0] < '0' || argv[2][0] > '9' || *end != '\0' ||
            errno != 0 || parsed == 0 || parsed > UINT32_MAX) {
            fputs("Seed must be a decimal integer from 1 to 4294967295.\n", stderr);
            return EXIT_FAILURE;
        }
        seed = (uint32_t)parsed;
    }
    if (!sollya_lib_init()) {
        fputs("Could not initialize libsollya\n", stderr);
        return EXIT_FAILURE;
    }
    FILE *output = fopen(output_path, "w");
    if (!output) {
        perror(output_path);
        sollya_lib_close();
        return EXIT_FAILURE;
    }
    FILE *details = fopen(detail_path, "w");
    if (!details) {
        perror(detail_path);
        fclose(output);
        sollya_lib_close();
        return EXIT_FAILURE;
    }

    fputs("degree\ttarget_condition\tactual_condition_median\t"
          "iter_comphorner_seconds\tsollya_seconds\tsuccessful_pairs\t"
          "iter_relative_error_lt_u_count\n", output);
    const char *detail_header =
        "degree\ttarget_condition\ttrial\tactual_condition\titer_relative_error\t"
        "sollya_relative_error\titer_status\titer_iterations\tsollya_status\t"
        "iter_faithful_verified\tsollya_faithful_verified\t"
        "iter_comphorner_seconds\tsollya_seconds\titer_relative_error_lt_u\t"
        "sollya_relative_error_lt_u\n";
    fputs(detail_header, stdout);
    fputs(detail_header, details);
    size_t cases = 0, successful = 0, unexpected = 0;
    size_t iter_below_u = 0, sollya_below_u = 0;
    uint32_t random_state = seed;
    fprintf(stderr, "Random seed: %" PRIu32 "\n", seed);
    double max_iter_error = 0.0, max_sollya_error = 0.0;
    for (size_t j = 0; j < sizeof(degrees) / sizeof(degrees[0]); ++j) {
        const size_t degree = degrees[j];
        for (int exponent = FIRST_EXPONENT; exponent <= LAST_EXPONENT; ++exponent) {
            double conditions[TRIALS], iter_times[TRIALS], sollya_times[TRIALS];
            size_t successful_pairs = 0;
            size_t iter_below_u_in_group = 0;
            bool complete = true;
            for (size_t trial = 0; trial < TRIALS; ++trial) {
                trial_result sample = {.condition = NAN, .iter_seconds = NAN,
                                       .sollya_seconds = NAN};
                const bool expected = run_trial(degree, exponent, trial, &random_state, &sample,
                                                &max_iter_error, &max_sollya_error, details);
                unexpected += !expected;
                successful_pairs += sample.both_faithful;
                iter_below_u_in_group += sample.iter_error_lt_u;
                sollya_below_u += sample.sollya_error_lt_u;
                conditions[trial] = sample.condition;
                iter_times[trial] = sample.iter_seconds;
                sollya_times[trial] = sample.sollya_seconds;
                complete &= isfinite(sample.condition) && isfinite(sample.iter_seconds) &&
                            isfinite(sample.sollya_seconds);
            }
            cases += TRIALS;
            successful += successful_pairs;
            iter_below_u += iter_below_u_in_group;
            fprintf(output, "%zu\t1e%d\t%.16e\t%.9e\t%.9e\t%zu\t%zu\n",
                    degree, exponent, complete ? median(conditions) : NAN,
                    complete ? median(iter_times) : NAN,
                    complete ? median(sollya_times) : NAN, successful_pairs, iter_below_u_in_group);
        }
    }

    sollya_lib_close();
    const bool summary_ok = close_report(output, output_path);
    const bool details_ok = close_report(details, detail_path);
    if (!summary_ok || !details_ok) {
        return EXIT_FAILURE;
    }
    fprintf(stderr, "%zu/%zu trials returned faithful results from both methods. "
            "Maximum relative error: "
            "IterCompHorner %.16e, Sollya (%d-bit) %.16e\n",
            successful, cases, max_iter_error, RESULT_PRECISION, max_sollya_error);
    fprintf(stderr, "Relative error < 2^-%d: IterCompHorner %zu/%zu, Sollya %zu/%zu\n",
            UNIT_ROUNDOFF_EXPONENT, iter_below_u, cases, sollya_below_u, cases);
    if (unexpected != 0) {
        fprintf(stderr, "%zu unexpected validation or timing failures\n", unexpected);
    }
    fprintf(stderr, "Wrote %s\nWrote %s\n", output_path, detail_path);
    return unexpected == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
