#include "functions.h"

#include <float.h>
#include <stdint.h>
#include <string.h>

#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024
#error "ufp_bits requires IEEE 754 binary64 doubles."
#endif

_Static_assert(sizeof(double) == sizeof(uint64_t),
               "ufp_bits requires IEEE 754 binary64 doubles.");

double ufp_bits(double x)
{
    const uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
    const uint64_t fraction_mask = UINT64_C(0x000fffffffffffff);
    uint64_t bits, exponent, fraction, result;
    double y;

    /* memcpy avoids undefined behaviour from pointer type punning. */
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
