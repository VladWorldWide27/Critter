#ifndef COMPUTATION_H
#define COMPUTATION_H

/*
 * computation - Critter prediction unit (compute-bound behaviour).
 *
 * Predicts HVAC behaviour by fitting a first-order thermal (RC) model to the
 * readings and extrapolating it. The model is Newton's law of cooling:
 *
 *     dT/dt = -k * (T - T_inf)
 *
 * where k is the thermal rate constant (1/s) and T_inf is the steady-state
 * temperature the room is heading toward. The fit is a nonlinear
 * least-squares problem, solved with Levenberg-Marquardt: each iteration
 * integrates the model with RK4, forms and solves a small normal-equation
 * system, and repeats to convergence. A Monte-Carlo ensemble of RK4
 * forecasts then yields an uncertainty band. This is deliberately heavier
 * than a single matrix multiply -- it is the compute-bound workload.
 */

#include <stddef.h>

#include "io.h"   /* critter_reading_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    COMP_NORMAL   = 0,  /* within the safe operating range                */
    COMP_WARNING  = 1,  /* at/above the warn threshold, or trending to it  */
    COMP_CRITICAL = 2   /* at/above the critical threshold                 */
} comp_status_t;

/* Tunables. Pass NULL to accept the defaults shown. */
typedef struct {
    double warn_c;      /* warning threshold, C   (default 30.0)  */
    double crit_c;      /* critical threshold, C  (default 35.0)  */
    double horizon_s;   /* forecast horizon, s    (default 300.0) */
    int    max_iter;    /* max LM iterations      (default 30)    */
    int    ensemble;    /* Monte-Carlo forecasts  (default 200)   */
    int    substeps;    /* RK4 substeps/interval  (default 4)     */
} comp_cfg_t;

typedef struct {
    double        current_c;        /* model temperature at the last sample */
    double        slope_c_per_s;    /* dT/dt now, from the fitted model      */
    double        forecast_c;       /* temperature at now + horizon_s        */
    double        forecast_hi_c;    /* ~95th-percentile forecast (ensemble)  */
    double        horizon_s;        /* horizon actually used                 */
    double        seconds_to_crit;  /* est. time to crit_c; < 0 if n/a       */
    int           will_breach;      /* 1 if forecast_c >= crit_c             */
    comp_status_t status;

    /* Fitted-model diagnostics. */
    double        k;                /* thermal rate constant, 1/s            */
    double        t_inf_c;          /* steady-state temperature, C           */
    double        rmse_c;           /* fit residual RMSE, C                  */
    int           iterations;       /* LM iterations actually run            */
} comp_result_t;

/* Analyse 'n' chronological readings. Returns 0 on success, -1 on invalid
 * arguments. With fewer than 3 readings the model is not identifiable and a
 * simple linear fallback is used. Allocation-free.
 *
 * The Monte-Carlo stage uses rand(); seed it (srand) beforehand if you want
 * varied ensembles across runs. */
int comp_analyze(const critter_reading_t *readings, size_t n,
                 const comp_cfg_t *cfg,
                 comp_result_t *out);

#ifdef __cplusplus
}
#endif

#endif /* COMPUTATION_H */