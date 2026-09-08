/*
 * mem - Critter memory-focused unit implementation.
 *
 * Portable C99. Link with -lm (sqrt, fabs, NAN).
 */

#include "mem.h"

#include <math.h>
#include <stdlib.h>

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a;
    double y = *(const double *)b;
    if (x < y) return -1;
    if (x > y) return  1;
    return 0;
}

static double median_sorted(const double *s, size_t n)
{
    if (n == 0)
        return 0.0;
    if (n & 1u)
        return s[n / 2];
    return 0.5 * (s[n / 2 - 1] + s[n / 2]);
}

/* Welford accumulator over the readings, optionally rejecting outliers. */
static void accumulate(const critter_reading_t *r, size_t n,
                       int reject, double median, double thresh,
                       double *mean_out, double *m2_out,
                       double *min_out, double *max_out,
                       size_t *kept_out, size_t *rej_out)
{
    double mean = 0.0, m2 = 0.0, mn = 0.0, mx = 0.0;
    size_t k = 0, rejected = 0, i;

    for (i = 0; i < n; i++) {
        double x = r[i].celsius;
        double d;

        if (reject && fabs(x - median) > thresh) {
            rejected++;
            continue;
        }
        k++;
        d     = x - mean;
        mean += d / (double)k;
        m2   += d * (x - mean);
        if (k == 1) {
            mn = mx = x;
        } else {
            if (x < mn) mn = x;
            if (x > mx) mx = x;
        }
    }

    *mean_out = mean;
    *m2_out   = m2;
    *min_out  = mn;
    *max_out  = mx;
    *kept_out = k;
    *rej_out  = rejected;
}

int mem_build(mem_summary_t *out,
              const critter_reading_t *r, size_t n,
              double *scratch,
              const mem_cfg_t *cfg)
{
    mem_cfg_t def = { 1, 3.5 };
    double median, mad, thresh, mean, m2, mn, mx;
    size_t k, rejected, i;
    int reject;

    if (!out || !r || !scratch || n == 0)
        return -1;
    if (!cfg)
        cfg = &def;

    /* median of the raw values */
    for (i = 0; i < n; i++)
        scratch[i] = r[i].celsius;
    qsort(scratch, n, sizeof(double), cmp_double);
    median = median_sorted(scratch, n);

    /* MAD = median of |x - median| (reuse scratch; order is irrelevant) */
    for (i = 0; i < n; i++)
        scratch[i] = fabs(scratch[i] - median);
    qsort(scratch, n, sizeof(double), cmp_double);
    mad = median_sorted(scratch, n);

    reject = cfg->reject_outliers && mad > 1e-12;
    thresh = cfg->mad_k * mad;

    accumulate(r, n, reject, median, thresh,
               &mean, &m2, &mn, &mx, &k, &rejected);

    /* Degenerate guard: if the threshold rejected everything, keep all. */
    if (k == 0) {
        accumulate(r, n, 0, median, thresh,
                   &mean, &m2, &mn, &mx, &k, &rejected);
        rejected = 0;
    }

    out->count_total    = n;
    out->count_inliers  = k;
    out->count_outliers = rejected;
    out->min_c          = mn;
    out->max_c          = mx;
    out->mean_c         = mean;
    out->variance_c     = k ? m2 / (double)k : 0.0;
    out->stddev_c       = sqrt(out->variance_c);
    out->median_c       = median;
    out->mad_c          = mad;
    out->t_start_ns     = r[0].timestamp_ns;
    out->t_end_ns       = r[n - 1].timestamp_ns;
    return 0;
}

size_t mem_decimate(const critter_reading_t *r, size_t n,
                    double *bucket_means, size_t nbuckets)
{
    uint64_t t0, t1, span;
    size_t b, i;

    if (!r || !bucket_means || n == 0 || nbuckets == 0)
        return 0;

    t0   = r[0].timestamp_ns;
    t1   = r[n - 1].timestamp_ns;
    span = (t1 > t0) ? (t1 - t0) : 0;

    for (b = 0; b < nbuckets; b++) {
        double sum = 0.0;
        size_t cnt = 0;

        for (i = 0; i < n; i++) {
            size_t bi;

            if (span == 0) {
                bi = 0;
            } else {
                bi = (size_t)(((r[i].timestamp_ns - t0) * (uint64_t)nbuckets)
                              / span);
                if (bi >= nbuckets)
                    bi = nbuckets - 1;
            }

            if (bi == b) {
                sum += r[i].celsius;
                cnt++;
            }
        }
        bucket_means[b] = cnt ? (sum / (double)cnt) : NAN;
    }
    return nbuckets;
}