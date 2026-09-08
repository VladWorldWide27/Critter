#ifndef IO_H
#define IO_H

/*
 * io - Critter I/O unit (I/O-bound behaviour).
 *
 * High-rate temperature acquisition into a fixed-capacity ring buffer.
 * No dynamic allocation: the caller supplies the storage. When the buffer
 * fills, the oldest reading is overwritten so the most recent window is
 * always retained for the downstream units.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* A single temperature reading. Produced by the I/O unit and consumed by the
 * mem (summary) and computation (prediction) units. */
typedef struct {
    double   celsius;       /* temperature in degrees Celsius             */
    uint64_t timestamp_ns;  /* CLOCK_MONOTONIC time of the sample, in ns  */
} critter_reading_t;

/* Pluggable temperature source. Return 0 on success and write the value to
 * *out_celsius; return non-zero on failure (the sample is then skipped). */
typedef int (*io_source_fn)(void *ctx, double *out_celsius);

/* High-rate acquisition unit. */
typedef struct {
    critter_reading_t *buf;         /* caller-owned, 'capacity' entries    */
    size_t             capacity;
    size_t             head;        /* next write index                    */
    size_t             count;       /* valid entries (<= capacity)         */
    uint64_t           total;       /* successful samples ever taken       */
    uint64_t           dropped;     /* source failures                     */
    io_source_fn       source;
    void              *source_ctx;
} io_t;

/* Initialise the unit over caller-provided storage. 'storage' must hold at
 * least 'capacity' readings and outlive the unit. If 'source' is NULL,
 * io_source_sysfs is used. Returns 0 on success, -1 on bad arguments. */
int io_init(io_t *io,
            critter_reading_t *storage, size_t capacity,
            io_source_fn source, void *source_ctx);

/* Take one sample now: read the source, timestamp it, push into the ring.
 * Returns 0 on success, -1 if the source failed. */
int io_sample_once(io_t *io);

/* Sample 'n_samples' times, paced at 'rate_hz' (samples/second) using an
 * absolute monotonic deadline so timing does not drift. This is the
 * I/O-bound work loop. rate_hz <= 0 means "as fast as possible".
 * Returns the number of successful samples. */
size_t io_run(io_t *io, double rate_hz, size_t n_samples);

/* Copy up to 'max' readings, oldest first, into 'dst' for the next stage.
 * Non-destructive. Returns the number copied. */
size_t io_snapshot(const io_t *io, critter_reading_t *dst, size_t max);

/* Number of readings currently buffered. */
size_t io_count(const io_t *io);

/* Drop all buffered readings (running counters are preserved). */
void io_reset(io_t *io);

/* Default source: reads a Linux thermal zone
 * (/sys/class/thermal/thermal_zone<n>/temp, millidegrees C).
 * 'ctx' may point to an int zone index, or be NULL for zone 0.
 * Present on Raspberry Pi 5 and most Linux hosts. */
int io_source_sysfs(void *ctx, double *out_celsius);

/* Convenience simulated source (slow drift + noise) for bring-up on hosts
 * without a usable thermal zone. 'ctx' is ignored. Not reentrant. */
int io_source_sim(void *ctx, double *out_celsius);

#ifdef __cplusplus
}
#endif

#endif /* IO_H */