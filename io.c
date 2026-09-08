/*
 * io - Critter I/O unit implementation.
 *
 * Portable C99. Uses POSIX clock_gettime / clock_nanosleep, available on
 * Raspberry Pi OS and on WSL/glibc. Link with -lm (simulated source uses sin).
 */

#define _POSIX_C_SOURCE 200809L

#include "io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
 * Timing back-end. glibc on Raspberry Pi OS and WSL provides CLOCK_MONOTONIC
 * and clock_nanosleep() once _POSIX_C_SOURCE (above) is set. If a toolchain --
 * or an IDE's code parser reading Windows headers -- does not expose them,
 * fall back to C-standard timing so the unit still builds. Define
 * IO_FORCE_PORTABLE to exercise the fallback deliberately.
 */
#if !defined(IO_FORCE_PORTABLE) && defined(CLOCK_MONOTONIC)
#  define IO_HAVE_POSIX_CLOCK 1
#else
#  define IO_HAVE_POSIX_CLOCK 0
#endif

static uint64_t now_ns(void)
{
#if IO_HAVE_POSIX_CLOCK
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
#else
    /* Coarser fallback (CPU-time based); adequate only where the POSIX
     * monotonic clock is unavailable. */
    return (uint64_t)((double)clock() * (1e9 / (double)CLOCKS_PER_SEC));
#endif
}

/* Sleep (or busy-wait) until the monotonic time reaches 'target_ns'. */
static void sleep_until_ns(uint64_t target_ns)
{
#if IO_HAVE_POSIX_CLOCK
    struct timespec d;
    d.tv_sec  = (time_t)(target_ns / 1000000000ull);
    d.tv_nsec = (long)(target_ns % 1000000000ull);
    clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &d, NULL);
#else
    while (now_ns() < target_ns) {
        /* busy-wait: ISO C has no portable sub-second sleep */
    }
#endif
}

int io_init(io_t *io,
            critter_reading_t *storage, size_t capacity,
            io_source_fn source, void *source_ctx)
{
    if (!io || !storage || capacity == 0)
        return -1;

    io->buf        = storage;
    io->capacity   = capacity;
    io->head       = 0;
    io->count      = 0;
    io->total      = 0;
    io->dropped    = 0;
    io->source     = source ? source : io_source_sysfs;
    io->source_ctx = source_ctx;
    return 0;
}

int io_sample_once(io_t *io)
{
    double c;

    if (!io)
        return -1;

    if (io->source(io->source_ctx, &c) != 0) {
        io->dropped++;
        return -1;
    }

    io->buf[io->head].celsius      = c;
    io->buf[io->head].timestamp_ns = now_ns();
    io->head = (io->head + 1) % io->capacity;

    if (io->count < io->capacity)
        io->count++;
    io->total++;
    return 0;
}

size_t io_run(io_t *io, double rate_hz, size_t n_samples)
{
    size_t ok = 0;
    size_t i;

    if (!io)
        return 0;

    if (rate_hz <= 0.0) {
        for (i = 0; i < n_samples; i++)
            if (io_sample_once(io) == 0)
                ok++;
        return ok;
    }

    {
        uint64_t period = (uint64_t)(1e9 / rate_hz);
        uint64_t next   = now_ns();

        for (i = 0; i < n_samples; i++) {
            if (io_sample_once(io) == 0)
                ok++;

            next += period;          /* absolute schedule -> no drift */
            sleep_until_ns(next);
        }
    }
    return ok;
}

size_t io_snapshot(const io_t *io, critter_reading_t *dst, size_t max)
{
    size_t n, start, i;

    if (!io || !dst)
        return 0;

    n     = io->count < max ? io->count : max;
    start = (io->head + io->capacity - io->count) % io->capacity;

    for (i = 0; i < n; i++)
        dst[i] = io->buf[(start + i) % io->capacity];
    return n;
}

size_t io_count(const io_t *io)
{
    return io ? io->count : 0;
}

void io_reset(io_t *io)
{
    if (io) {
        io->head  = 0;
        io->count = 0;
    }
}

int io_source_sysfs(void *ctx, double *out_celsius)
{
    int  zone = ctx ? *(int *)ctx : 0;
    char path[64];
    FILE *f;
    long milli;

    if (!out_celsius)
        return -1;

    snprintf(path, sizeof path,
             "/sys/class/thermal/thermal_zone%d/temp", zone);

    f = fopen(path, "r");
    if (!f)
        return -1;

    if (fscanf(f, "%ld", &milli) != 1) {
        fclose(f);
        return -1;
    }
    fclose(f);

    *out_celsius = (double)milli / 1000.0;
    return 0;
}

int io_source_sim(void *ctx, double *out_celsius)
{
    static double phase = 0.0;
    double noise;

    (void)ctx;
    if (!out_celsius)
        return -1;

    phase += 0.05;
    noise  = ((double)rand() / (double)RAND_MAX - 0.5) * 0.4;

    /* ~22 C baseline, slow +/-3 C drift, small pseudo-random noise */
    *out_celsius = 22.0 + 3.0 * sin(phase) + noise;
    return 0;
}