#ifndef MEM_H
#define MEM_H

/*
 * mem - Critter memory-focused unit (memory-bound behaviour).
 *
 * Reduces an arbitrarily long capture into a small, fixed-size form for
 * off-line analysis: robust statistics with outlier rejection, plus an
 * optional decimation of the series into a handful of time buckets.
 * No dynamic allocation: the caller supplies scratch/output storage.
 */

#include <stddef.h>
#include <stdint.h>

#include "io.h"   /* critter_reading_t */

#ifdef __cplusplus
extern "C" {
#endif

/* Compact description of a batch of readings. Independent of input size. */
typedef struct {
    size_t   count_total;    /* readings examined                          */
    size_t   count_inliers;  /* kept after outlier rejection               */
    size_t   count_outliers; /* rejected                                   */
    double   min_c;          /* stats below are over the inliers           */
    double   max_c;
    double   mean_c;
    double   variance_c;     /* population variance                        */
    double   stddev_c;
    double   median_c;       /* over all readings                          */
    double   mad_c;          /* median absolute deviation, raw scale       */
    uint64_t t_start_ns;     /* timestamp of first reading                 */
    uint64_t t_end_ns;       /* timestamp of last reading                  */
} mem_summary_t;

/* Tunables. Pass NULL to any API to accept the defaults. */
typedef struct {
    int    reject_outliers;  /* non-zero enables MAD rejection (default 1) */
    double mad_k;            /* reject |x - median| > k*MAD (default 3.5)  */
} mem_cfg_t;

/* Summarise 'n' readings with robust (MAD-based) outlier rejection.
 * 'scratch' must provide room for 'n' doubles and is used internally so the
 * unit allocates nothing. Returns 0 on success, -1 on invalid arguments. */
int mem_build(mem_summary_t *out,
              const critter_reading_t *readings, size_t n,
              double *scratch,
              const mem_cfg_t *cfg);

/* Reduce a series to 'nbuckets' equal-time buckets, writing each bucket's
 * mean to bucket_means[] (empty buckets become NaN). This is the memory
 * optimisation: a long capture collapses to a fixed footprint.
 * Cost is O(n * nbuckets); nbuckets is expected small.
 * Returns nbuckets on success, 0 on invalid arguments. */
size_t mem_decimate(const critter_reading_t *readings, size_t n,
                    double *bucket_means, size_t nbuckets);

#ifdef __cplusplus
}
#endif

#endif /* MEM_H */