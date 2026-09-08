#include "io.h"
#include "mem.h"
#include "computation.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    srand(1234);                               /* varies the MC ensemble */

    critter_reading_t ring[512], snap[512];
    io_t io;
    io_init(&io, ring, 512, io_source_sim, NULL);   /* sim source; use NULL for the Pi sensor */
    io_run(&io, 2.0, 120);                          /* 2 Hz for 120 samples = 60 s */

    size_t n = io_snapshot(&io, snap, 512);

    double scratch[512], buckets[8];
    mem_summary_t s;
    mem_build(&s, snap, n, scratch, NULL);
    mem_decimate(snap, n, buckets, 8);
    printf("mean=%.2fC  sd=%.2f  outliers=%zu/%zu\n",
           s.mean_c, s.stddev_c, s.count_outliers, s.count_total);

    comp_result_t p;
    comp_analyze(snap, n, NULL, &p);
    const char *st = p.status == COMP_CRITICAL ? "CRITICAL"
                   : p.status == COMP_WARNING  ? "WARNING" : "NORMAL";
    printf("now=%.2fC  fit k=%.4f Tinf=%.2f  forecast(+%.0fs)=%.2f  %s\n",
           p.current_c, p.k, p.t_inf_c, p.horizon_s, p.forecast_c, st);
    return 0;
}