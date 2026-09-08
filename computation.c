/*
 * computation - Critter prediction unit implementation.
 *
 * Nonlinear least-squares fit of a first-order thermal model via
 * Levenberg-Marquardt, using RK4 for the model response and a numeric
 * (finite-difference) Jacobian, followed by an RK4 Monte-Carlo forecast
 * ensemble. Two parameters are estimated: k and T_inf. The initial
 * temperature is anchored to the first reading.
 *
 * Portable C99. Link with -lm.
 */

#include "computation.h"

#include <math.h>
#include <stdlib.h>

/* ---- model + integrator ------------------------------------------------ */

/* Newton's law of cooling. */
static double f_ode(double T, double k, double t_inf)
{
    return -k * (T - t_inf);
}

/* Advance T over an interval of length h using RK4 with 'sub' substeps. */
static double rk4_advance(double T, double k, double t_inf, double h, int sub)
{
    double hh = h / (double)sub;
    int s;

    if (sub < 1) sub = 1;
    for (s = 0; s < sub; s++) {
        double k1 = f_ode(T,               k, t_inf);
        double k2 = f_ode(T + 0.5 * hh * k1, k, t_inf);
        double k3 = f_ode(T + 0.5 * hh * k2, k, t_inf);
        double k4 = f_ode(T +       hh * k3, k, t_inf);
        T += (hh / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    return T;
}

/* Sum of squared residuals for a parameter set (single base trajectory). */
static double model_cost(double k, double t_inf, double T0,
                         const critter_reading_t *r, size_t n, int sub)
{
    double T = T0;
    double prev_t = (double)r[0].timestamp_ns / 1e9;
    double cost = 0.0;
    size_t i;

    /* residual at i = 0 is zero: the model is anchored to r[0]. */
    for (i = 1; i < n; i++) {
        double t = (double)r[i].timestamp_ns / 1e9;
        double h = t - prev_t;
        double d;
        prev_t = t;
        if (h < 0.0) h = 0.0;
        T = rk4_advance(T, k, t_inf, h, sub);
        d = T - r[i].celsius;
        cost += d * d;
    }
    return cost;
}

/* Build the 2x2 normal-equation system J^T J and J^T r for the current
 * parameters, using forward-difference sensitivities. The base and both
 * perturbed trajectories are integrated in lockstep, so no O(n) storage is
 * needed. */
static void model_normal_eqs(const double theta[2], double T0,
                             const critter_reading_t *r, size_t n, int sub,
                             double JtJ[2][2], double Jtr[2], double *cost_out)
{
    double Tb  = T0;              /* base trajectory                      */
    double Tp0 = T0, Tp1 = T0;    /* perturbed in k, perturbed in t_inf   */
    double hk  = (fabs(theta[0]) > 1e-8) ? 1e-4 * fabs(theta[0]) : 1e-6;
    double hi  = (fabs(theta[1]) > 1e-8) ? 1e-4 * fabs(theta[1]) : 1e-6;
    double prev_t = (double)r[0].timestamp_ns / 1e9;
    double cost = 0.0;
    size_t i;

    JtJ[0][0] = JtJ[0][1] = JtJ[1][0] = JtJ[1][1] = 0.0;
    Jtr[0] = Jtr[1] = 0.0;

    for (i = 1; i < n; i++) {
        double t = (double)r[i].timestamp_ns / 1e9;
        double h = t - prev_t;
        double res, d0, d1;
        prev_t = t;
        if (h < 0.0) h = 0.0;

        Tb  = rk4_advance(Tb,  theta[0],      theta[1],      h, sub);
        Tp0 = rk4_advance(Tp0, theta[0] + hk, theta[1],      h, sub);
        Tp1 = rk4_advance(Tp1, theta[0],      theta[1] + hi, h, sub);

        res = Tb - r[i].celsius;
        cost += res * res;

        d0 = (Tp0 - Tb) / hk;    /* d(pred)/dk     */
        d1 = (Tp1 - Tb) / hi;    /* d(pred)/dt_inf */

        JtJ[0][0] += d0 * d0;
        JtJ[0][1] += d0 * d1;
        JtJ[1][1] += d1 * d1;
        Jtr[0]    += d0 * res;
        Jtr[1]    += d1 * res;
    }
    JtJ[1][0] = JtJ[0][1];
    *cost_out = cost;
}

/* Invert a 2x2 matrix. Returns 0 on success, -1 if (near-)singular. */
static int inv2(double A[2][2], double out[2][2])
{
    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (fabs(det) < 1e-300)
        return -1;
    {
        double inv = 1.0 / det;
        out[0][0] =  A[1][1] * inv;
        out[0][1] = -A[0][1] * inv;
        out[1][0] = -A[1][0] * inv;
        out[1][1] =  A[0][0] * inv;
    }
    return 0;
}

/* Standard-normal sample via Box-Muller, one cached value between calls. */
static double gauss_rand(void)
{
    static int    have = 0;
    static double cached = 0.0;
    double u1, u2, mag;

    if (have) {
        have = 0;
        return cached;
    }
    do {
        u1 = (double)rand() / (double)RAND_MAX;
    } while (u1 <= 1e-12);
    u2 = (double)rand() / (double)RAND_MAX;

    mag    = sqrt(-2.0 * log(u1));
    cached = mag * sin(2.0 * 3.14159265358979323846 * u2);
    have   = 1;
    return mag * cos(2.0 * 3.14159265358979323846 * u2);
}

/* ---- public API -------------------------------------------------------- */

int comp_analyze(const critter_reading_t *r, size_t n,
                 const comp_cfg_t *cfg, comp_result_t *out)
{
    comp_cfg_t def = { 30.0, 35.0, 300.0, 30, 200, 4 };
    double T0, span, theta[2], cost, lambda, mean;
    double current, slope, forecast, forecast_hi, tt_crit;
    double sigma, sd_k, sd_i;
    int sub, it, iters, breach;
    size_t i;
    comp_status_t status;

    if (!out || !r || n == 0)
        return -1;
    if (!cfg)
        cfg = &def;

    sub  = cfg->substeps > 0 ? cfg->substeps : 4;
    T0   = r[0].celsius;
    span = (double)r[n - 1].timestamp_ns / 1e9 - (double)r[0].timestamp_ns / 1e9;
    if (span < 0.0) span = 0.0;

    /* ---- fewer than 3 points: model not identifiable, linear fallback -- */
    if (n < 3) {
        current = r[n - 1].celsius;
        slope   = 0.0;
        if (n == 2 && span > 1e-9)
            slope = (r[1].celsius - r[0].celsius) / span;

        forecast    = current + slope * cfg->horizon_s;
        forecast_hi = forecast;
        tt_crit = (slope > 1e-9 && current < cfg->crit_c)
                ? (cfg->crit_c - current) / slope : -1.0;

        out->k          = 0.0;
        out->t_inf_c    = current;
        out->rmse_c     = 0.0;
        out->iterations = 0;
        goto finish;
    }

    /* ---- initial guess ------------------------------------------------- */
    mean = 0.0;
    for (i = 0; i < n; i++)
        mean += r[i].celsius;
    mean /= (double)n;

    theta[0] = 0.01;   /* k     */
    theta[1] = mean;   /* t_inf */

    /* ---- Levenberg-Marquardt ------------------------------------------ */
    cost   = model_cost(theta[0], theta[1], T0, r, n, sub);
    lambda = 1e-2;
    iters  = 0;

    for (it = 0; it < cfg->max_iter; it++) {
        double JtJ[2][2], Jtr[2], c_here;
        int accepted = 0, tries;

        model_normal_eqs(theta, T0, r, n, sub, JtJ, Jtr, &c_here);

        for (tries = 0; tries < 8; tries++) {
            double A[2][2], inv[2][2], step[2], cand[2], cnew;

            A[0][0] = JtJ[0][0] * (1.0 + lambda) + 1e-9;
            A[1][1] = JtJ[1][1] * (1.0 + lambda) + 1e-9;
            A[0][1] = JtJ[0][1];
            A[1][0] = JtJ[1][0];

            if (inv2(A, inv) != 0) {
                lambda *= 4.0;
                continue;
            }

            step[0] = -(inv[0][0] * Jtr[0] + inv[0][1] * Jtr[1]);
            step[1] = -(inv[1][0] * Jtr[0] + inv[1][1] * Jtr[1]);

            cand[0] = theta[0] + step[0];
            cand[1] = theta[1] + step[1];
            if (cand[0] < 1e-6) cand[0] = 1e-6;   /* keep k > 0 */

            cnew = model_cost(cand[0], cand[1], T0, r, n, sub);

            if (cnew < cost) {
                double rel = (cost - cnew) / (cost > 1e-30 ? cost : 1e-30);
                theta[0] = cand[0];
                theta[1] = cand[1];
                cost     = cnew;
                lambda  *= 0.5;
                if (lambda < 1e-9) lambda = 1e-9;
                iters++;
                accepted = (rel < 1e-8) ? 2 : 1;   /* 2 == converged */
                break;
            }
            lambda *= 4.0;
            if (lambda > 1e12)
                break;
        }

        if (accepted == 0)   /* no downhill step found */
            break;
        if (accepted == 2)   /* relative improvement below tolerance */
            break;
    }

    /* ---- point prediction from the fitted model ------------------------ */
    current  = theta[1] + (T0 - theta[1]) * exp(-theta[0] * span);
    slope    = -theta[0] * (current - theta[1]);
    forecast = theta[1] + (current - theta[1]) * exp(-theta[0] * cfg->horizon_s);

    tt_crit = -1.0;
    {
        double num = cfg->crit_c - theta[1];
        double den = current      - theta[1];
        if (theta[0] > 1e-9 && fabs(den) > 1e-12) {
            double ratio = num / den;
            if (ratio > 0.0 && ratio < 1.0) {
                double tau = -log(ratio) / theta[0];
                if (isfinite(tau) && tau >= 0.0)
                    tt_crit = tau;
            }
        }
    }

    /* ---- Monte-Carlo forecast ensemble --------------------------------- */
    forecast_hi = forecast;
    {
        double JtJ[2][2], Jtr[2], c_final, invJ[2][2];
        int M = cfg->ensemble;

        model_normal_eqs(theta, T0, r, n, sub, JtJ, Jtr, &c_final);
        sigma = sqrt(c_final / (double)((n > 2) ? (n - 2) : 1));

        sd_k = sd_i = 0.0;
        if (inv2(JtJ, invJ) == 0) {
            double v0 = sigma * sigma * invJ[0][0];
            double v1 = sigma * sigma * invJ[1][1];
            if (v0 > 0.0) sd_k = sqrt(v0);
            if (v1 > 0.0) sd_i = sqrt(v1);
        }

        if (M > 0 && (sd_k > 0.0 || sd_i > 0.0)) {
            double fmean = 0.0, fm2 = 0.0;
            int m;
            for (m = 0; m < M; m++) {
                double kk = theta[0] + sd_k * gauss_rand();
                double ti = theta[1] + sd_i * gauss_rand();
                double fc, d;
                if (kk < 1e-6) kk = 1e-6;
                fc = rk4_advance(current, kk, ti, cfg->horizon_s, 16);
                d  = fc - fmean;
                fmean += d / (double)(m + 1);
                fm2   += d * (fc - fmean);
            }
            {
                double fsd = sqrt(fm2 / (double)M);
                forecast_hi = fmean + 1.645 * fsd;   /* ~95th percentile */
            }
        }
    }

    out->k          = theta[0];
    out->t_inf_c    = theta[1];
    out->rmse_c     = sqrt(cost / (double)n);
    out->iterations = iters;

finish:
    breach = (forecast >= cfg->crit_c) ? 1 : 0;

    if (current >= cfg->crit_c)
        status = COMP_CRITICAL;
    else if (current >= cfg->warn_c)
        status = COMP_WARNING;
    else if (breach || forecast_hi >= cfg->crit_c || forecast >= cfg->warn_c)
        status = COMP_WARNING;
    else
        status = COMP_NORMAL;

    out->current_c       = current;
    out->slope_c_per_s   = slope;
    out->forecast_c      = forecast;
    out->forecast_hi_c   = forecast_hi;
    out->horizon_s       = cfg->horizon_s;
    out->seconds_to_crit = tt_crit;
    out->will_breach     = breach;
    out->status          = status;
    return 0;
}