#include "ewcr.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MMC_MAX_THETA 12
#define MMC_MAX_MODELS 40
#define MMC_MAX_N 512

typedef struct {
    char name[40];
    char family[20];
    char base_family[20];
    int order, eta, factor, kx, p;
    double theta_hat[MMC_MAX_THETA];
    double theta_center[MMC_MAX_THETA];
    double theta_width[MMC_MAX_THETA];
} MmcModel;

typedef struct {
    char name[8];
    char transform[8];
    int nbits, kr, value_bits, nsel;
    int indices[MMC_MAX_N];
    int q[MMC_MAX_N];
    double rhat[MMC_MAX_N];
} MmcResidual;

typedef struct {
    MmcModel model;
    MmcResidual residual;
    int kx, nx, nr, header_bits, total_bits, nm, nl;
    double theta_q[MMC_MAX_THETA];
    double mse;
    double xhat[MMC_MAX_N];
} MmcPacket;

typedef struct {
    double prev2[MMC_MAX_N];
    double prev1[MMC_MAX_N];
    int has_prev_model;
    MmcModel prev_model;
    double prev_theta_q[MMC_MAX_THETA];
    int history_count;
} MmcState;

static double poly_width_data[10][10];
static int poly_width_ready = 0;

static void init_poly_width(void)
{
    if (poly_width_ready) return;
    double sigma[10][10] = {
        {0.215},
        {0.2150, 0.3812},
        {0.2409, 0.3812, 0.4191},
        {0.2409, 0.2812, 0.4192, 0.2539},
        {0.2258, 0.2812, 0.3572, 0.2539, 0.1363},
        {0.2258, 0.2877, 0.3572, 0.2310, 0.1364, 0.0654},
        {0.2268, 0.2877, 0.3611, 0.2310, 0.1264, 0.0655, 0.0426},
        {0.2268, 0.2870, 0.3611, 0.2316, 0.1264, 0.0601, 0.0427, 0.0326},
        {0.2268, 0.2871, 0.3609, 0.2316, 0.1262, 0.0601, 0.0383, 0.0327, 0.0280},
        {0.2268, 0.2871, 0.3609, 0.2317, 0.1261, 0.0593, 0.0384, 0.0289, 0.0281, 0.0244}
    };
    int lens[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    for (int k = 0; k < 10; k++) {
        for (int j = 0; j < lens[k]; j++) {
            double w = 7.0 * sigma[k][j];
            if (w > 1.55) w = 1.55;
            poly_width_data[k][j] = w;
        }
    }
    poly_width_ready = 1;
}

static void cheby_basis(int N, int ord, double *B)
{
    /* B is N x (ord+1) row-major */
    int p = ord + 1;
    for (int i = 0; i < N; i++) {
        double t = -1.0 + (1.0 - 1.0 / N + 1.0) * i / (N - 1 == 0 ? 1 : (N - 1));
        /* linspace(-1, 1-1/N, N) */
        t = -1.0 + ((1.0 - 1.0 / N) - (-1.0)) * i / (double)(N - 1);
        B[i * p + 0] = 1.0;
        if (ord >= 1) B[i * p + 1] = t;
        for (int k = 2; k <= ord; k++)
            B[i * p + k] = 2.0 * t * B[i * p + (k - 1)] - B[i * p + (k - 2)];
    }
}

static void predictor_matrix(const double *xprev, int nprev, int N, int ord, int eta, double *X)
{
    memset(X, 0, (size_t)N * ord * sizeof(double));
    for (int i = 0; i < N; i++) {
        int base = N + (i + 1) - eta; /* MATLAB i=1..N */
        for (int j = 1; j <= ord; j++) {
            int idx = base - j + 1; /* 1-based */
            if (idx >= 1 && idx <= nprev) X[i * ord + (j - 1)] = xprev[idx - 1];
        }
    }
}

static void reconstruct_model(const MmcModel *m, const double *theta, const MmcState *st, int N,
                              double fs, double *y)
{
    const char *fam = m->family;
    if (strcmp(fam, "pred_para") == 0 && m->base_family[0]) fam = m->base_family;
    memset(y, 0, (size_t)N * sizeof(double));
    if (strcmp(fam, "none") == 0) return;
    if (strcmp(fam, "sin") == 0) {
        for (int i = 0; i < N; i++) {
            double t = i / fs;
            y[i] = theta[0] * cos(2 * M_PI * theta[1] * t + theta[2]);
        }
        return;
    }
    if (strcmp(fam, "poly") == 0) {
        int p = m->order + 1;
        double *B = (double *)ewcr_xmalloc((size_t)N * p * sizeof(double));
        cheby_basis(N, m->order, B);
        for (int i = 0; i < N; i++) {
            double s = 0.0;
            for (int j = 0; j < p; j++) s += B[i * p + j] * theta[j];
            y[i] = s;
        }
        free(B);
        return;
    }
    if (strcmp(fam, "pred_samples") == 0) {
        double xprev[2 * MMC_MAX_N];
        double scale = pow(2.0, -m->kx);
        for (int i = 0; i < N; i++) {
            xprev[i] = st->prev2[i] * scale;
            xprev[N + i] = st->prev1[i] * scale;
        }
        int ord = m->order;
        if (ord < 1) return;
        double *X = (double *)ewcr_xmalloc((size_t)N * ord * sizeof(double));
        predictor_matrix(xprev, 2 * N, N, ord, m->eta, X);
        for (int i = 0; i < N; i++) {
            double s = 0.0;
            for (int j = 0; j < ord; j++) s += X[i * ord + j] * theta[j];
            y[i] = s;
        }
        free(X);
    }
}

static void quantize_vector(const double *theta, const int *bits, const double *center,
                            const double *width, int p, double *q, int *idx)
{
    for (int j = 0; j < p; j++) {
        int b = bits[j];
        double c = center[j], w = width[j];
        if (b <= 0 || w <= 0) {
            idx[j] = 0;
            q[j] = c;
        } else {
            double delta = w / (double)(1 << b);
            int ind = (int)floor((theta[j] - c) / delta);
            int him = (1 << (b - 1)) - 1;
            int lom = -(1 << (b - 1));
            if (ind > him) ind = him;
            if (ind < lom) ind = lom;
            idx[j] = ind;
            q[j] = delta * (ind + 0.5) + c;
        }
    }
}

static void allocate_bits_greedy(const MmcModel *m, int nx, const MmcState *st, int N, double fs,
                                 int minb, int maxb, int *bits)
{
    int p = m->p;
    for (int j = 0; j < p; j++) bits[j] = minb;
    int sum = p * minb;
    if (sum > nx) {
        for (int j = 0; j < p; j++) bits[j] = 0;
        sum = 0;
    }
    double yun[MMC_MAX_N];
    reconstruct_model(m, m->theta_hat, st, N, fs, yun);
    while (sum < nx) {
        int allmax = 1;
        for (int j = 0; j < p; j++) if (bits[j] < maxb) allmax = 0;
        if (allmax) break;
        double q0[MMC_MAX_THETA];
        int idx0[MMC_MAX_THETA];
        quantize_vector(m->theta_hat, bits, m->theta_center, m->theta_width, p, q0, idx0);
        double y0[MMC_MAX_N];
        reconstruct_model(m, q0, st, N, fs, y0);
        double e0 = 0.0;
        for (int i = 0; i < N; i++) {
            double d = yun[i] - y0[i];
            e0 += d * d;
        }
        double bestGain = -1e300;
        int bestj = 0;
        for (int j = 0; j < p; j++) {
            if (bits[j] >= maxb) continue;
            int bt[MMC_MAX_THETA];
            memcpy(bt, bits, (size_t)p * sizeof(int));
            bt[j]++;
            double q[MMC_MAX_THETA];
            int id[MMC_MAX_THETA];
            quantize_vector(m->theta_hat, bt, m->theta_center, m->theta_width, p, q, id);
            double yh[MMC_MAX_N];
            reconstruct_model(m, q, st, N, fs, yh);
            double e = 0.0;
            for (int i = 0; i < N; i++) {
                double d = yun[i] - yh[i];
                e += d * d;
            }
            double gain = e0 - e;
            if (gain > bestGain) {
                bestGain = gain;
                bestj = j;
            }
        }
        bits[bestj]++;
        sum++;
    }
}

static void fit_sinusoid(const double *y, int N, double fs, const double *center, const double *width,
                         double *theta)
{
    memcpy(theta, center, 3 * sizeof(double));
    double t0 = (double)N;
    (void)t0;
    double flo = center[1] - width[1] / 2.0;
    double fhi = center[1] + width[1] / 2.0;
    int ngrid = (fhi <= flo) ? 1 : 81;
    double best = 1e300;
    double *A = (double *)ewcr_xmalloc((size_t)N * 2 * sizeof(double));
    double ab[2];
    for (int g = 0; g < ngrid; g++) {
        double f = (ngrid == 1) ? center[1] : flo + (fhi - flo) * g / 80.0;
        for (int i = 0; i < N; i++) {
            double t = i / fs;
            A[i * 2 + 0] = cos(2 * M_PI * f * t);
            A[i * 2 + 1] = sin(2 * M_PI * f * t);
        }
        ewcr_ls_solve(A, N, 2, y, ab);
        double amp = hypot(ab[0], ab[1]);
        double phi = atan2(-ab[1], ab[0]);
        double th[3] = {amp, f, phi};
        for (int k = 0; k < 3; k++) {
            double lo = center[k] - width[k] / 2.0;
            double hi = center[k] + width[k] / 2.0;
            if (th[k] < lo) th[k] = lo;
            if (th[k] > hi) th[k] = hi;
        }
        double e = 0.0;
        for (int i = 0; i < N; i++) {
            double t = i / fs;
            double yh = th[0] * cos(2 * M_PI * th[1] * t + th[2]);
            double d = y[i] - yh;
            e += d * d;
        }
        if (e < best) {
            best = e;
            memcpy(theta, th, 3 * sizeof(double));
        }
    }
    free(A);
}

static void fit_poly(const double *y, int N, int ord, const double *c, const double *w, double *theta)
{
    int p = ord + 1;
    double *B = (double *)ewcr_xmalloc((size_t)N * p * sizeof(double));
    cheby_basis(N, ord, B);
    ewcr_ls_solve(B, N, p, y, theta);
    for (int j = 0; j < p; j++) {
        double lo = c[j] - w[j] / 2.0, hi = c[j] + w[j] / 2.0;
        if (theta[j] < lo) theta[j] = lo;
        if (theta[j] > hi) theta[j] = hi;
    }
    free(B);
}

static void fit_sample_pred(const double *y, const double *xprev, int nprev, int N, int ord, int eta,
                            const double *c, const double *w, double *theta)
{
    double *X = (double *)ewcr_xmalloc((size_t)N * ord * sizeof(double));
    predictor_matrix(xprev, nprev, N, ord, eta, X);
    double fro = 0.0, ny = 0.0;
    for (int i = 0; i < N * ord; i++) fro += X[i] * X[i];
    for (int i = 0; i < N; i++) ny += y[i] * y[i];
    fro = sqrt(fro);
    if (fro <= 1e-14 * fmax(1.0, sqrt(ny))) {
        memset(theta, 0, (size_t)ord * sizeof(double));
    } else {
        ewcr_ls_solve(X, N, ord, y, theta);
    }
    for (int j = 0; j < ord; j++) {
        double lo = c[j] - w[j] / 2.0, hi = c[j] + w[j] / 2.0;
        if (theta[j] < lo) theta[j] = lo;
        if (theta[j] > hi) theta[j] = hi;
    }
    free(X);
}

static MmcModel mk(const char *name, const char *family, int order, int eta, int factor,
                   const double *theta, int p, const double *c, const double *w, int kx)
{
    MmcModel s;
    memset(&s, 0, sizeof(s));
    strncpy(s.name, name, sizeof(s.name) - 1);
    strncpy(s.family, family, sizeof(s.family) - 1);
    s.order = order;
    s.eta = eta;
    s.factor = factor;
    s.kx = kx;
    s.p = p;
    for (int i = 0; i < p; i++) {
        s.theta_hat[i] = theta ? theta[i] : 0;
        s.theta_center[i] = c ? c[i] : 0;
        s.theta_width[i] = w ? w[i] : 0;
    }
    return s;
}

static void encode_residual(const double *r, int N, const char *method, int budget, int n_kr,
                            int value_bits, int min_coeffs, MmcResidual *res)
{
    memset(res, 0, sizeof(*res));
    strncpy(res->name, "none", sizeof(res->name) - 1);
    res->value_bits = value_bits;
    if (strcmp(method, "none") == 0 || strcmp(method, "NONE") == 0 || budget <= 0) return;

    double *c = (double *)ewcr_xmalloc((size_t)N * sizeof(double));
    if (strcmp(method, "DCT") == 0) ewcr_dct_fwd(r, c, N);
    else if (strcmp(method, "DWT") == 0) ewcr_haar_fwd(r, c, N);
    else {
        free(c);
        return;
    }
    double mx = 0.0;
    for (int i = 0; i < N; i++) {
        double a = fabs(c[i]);
        if (a > mx) mx = a;
    }
    int kr = (int)ceil(log2(mx + 1e-12));
    if (kr < 0) kr = 0;
    int krmax = (1 << n_kr) - 1;
    if (kr > krmax) kr = krmax;
    double *cn = (double *)ewcr_xmalloc((size_t)N * sizeof(double));
    double sc = pow(2.0, -kr);
    for (int i = 0; i < N; i++) cn[i] = c[i] * sc;

    int index_bits = (int)ceil(log2((double)N));
    int bits_per = index_bits + value_bits;
    int K = budget / bits_per;
    if (K > N) K = N;
    if (K < 0) K = 0;
    strncpy(res->name, method, sizeof(res->name) - 1);
    strncpy(res->transform, method, sizeof(res->transform) - 1);
    res->kr = kr;
    if (K < min_coeffs) {
        free(c);
        free(cn);
        return;
    }
    int *ord = (int *)ewcr_xmalloc((size_t)N * sizeof(int));
    for (int i = 0; i < N; i++) ord[i] = i;
    for (int i = 0; i < K; i++) { /* partial sort top K */
        int bi = i;
        for (int j = i + 1; j < N; j++) if (fabs(cn[ord[j]]) > fabs(cn[ord[bi]])) bi = j;
        int t = ord[i];
        ord[i] = ord[bi];
        ord[bi] = t;
    }
    int *sel = (int *)ewcr_xmalloc((size_t)K * sizeof(int));
    for (int i = 0; i < K; i++) sel[i] = ord[i];
    for (int i = 1; i < K; i++) {
        int v = sel[i], j = i - 1;
        while (j >= 0 && sel[j] > v) {
            sel[j + 1] = sel[j];
            j--;
        }
        sel[j + 1] = v;
    }
    int qmax = (1 << (value_bits - 1)) - 1;
    double *chat = (double *)ewcr_xcalloc((size_t)N, sizeof(double));
    res->nsel = K;
    for (int i = 0; i < K; i++) {
        int qi = (int)round(cn[sel[i]] * qmax);
        if (qi > qmax) qi = qmax;
        if (qi < -qmax) qi = -qmax;
        res->indices[i] = sel[i];
        res->q[i] = qi;
        chat[sel[i]] = (double)qi / qmax;
    }
    double *c_hat = (double *)ewcr_xmalloc((size_t)N * sizeof(double));
    double p2 = pow(2.0, kr);
    for (int i = 0; i < N; i++) c_hat[i] = chat[i] * p2;
    if (strcmp(method, "DCT") == 0) ewcr_dct_inv(c_hat, res->rhat, N);
    else ewcr_haar_inv(c_hat, res->rhat, N);
    res->nbits = K * bits_per;
    free(c);
    free(cn);
    free(ord);
    free(sel);
    free(chat);
    free(c_hat);
}

static void decode_residual(const MmcResidual *res, int N, double *rhat)
{
    memset(rhat, 0, (size_t)N * sizeof(double));
    if (strcmp(res->name, "none") == 0 || res->nbits == 0) return;
    int qmax = (1 << (res->value_bits - 1)) - 1;
    double *c = (double *)ewcr_xcalloc((size_t)N, sizeof(double));
    double p2 = pow(2.0, res->kr);
    for (int i = 0; i < res->nsel; i++) {
        int idx = res->indices[i];
        if (idx >= 0 && idx < N) c[idx] = ((double)res->q[i] / qmax) * p2;
    }
    if (strcmp(res->transform, "DCT") == 0) ewcr_dct_inv(c, rhat, N);
    else if (strcmp(res->transform, "DWT") == 0) ewcr_haar_inv(c, rhat, N);
    free(c);
}

static int encode_window(const double *x, int N, MmcState *st, const MmcOpts *cfg, MmcPacket *best)
{
    init_poly_width();
    double mx = 0.0;
    for (int i = 0; i < N; i++) {
        double a = fabs(x[i]);
        if (a > mx) mx = a;
    }
    int kx = (int)ceil(log2(mx + 1e-8));
    if (kx < 0) kx = 0;
    int kxmax = (1 << cfg->n_kx) - 1;
    if (kx > kxmax) kx = kxmax;
    double *xn = (double *)ewcr_xmalloc((size_t)N * sizeof(double));
    double sc = pow(2.0, -kx);
    for (int i = 0; i < N; i++) xn[i] = x[i] * sc;

    MmcModel models[MMC_MAX_MODELS];
    int nm = 0;
    models[nm++] = mk("none", "none", 0, 0, 0, NULL, 0, NULL, NULL, kx);

    {
        double c[3] = {0.75, cfg->fn, 0};
        double w[3] = {0.5, 0.2, 2 * M_PI};
        double th[3];
        fit_sinusoid(xn, N, cfg->fs, c, w, th);
        models[nm++] = mk("sin-1", "sin", 0, 0, 0, th, 3, c, w, kx);
    }
    for (int ord = 0; ord <= 9; ord++) {
        double c[MMC_MAX_THETA], w[MMC_MAX_THETA], th[MMC_MAX_THETA];
        int p = ord + 1;
        for (int j = 0; j < p; j++) {
            c[j] = 0;
            w[j] = poly_width_data[ord][j];
        }
        fit_poly(xn, N, ord, c, w, th);
        char name[32];
        snprintf(name, sizeof(name), "poly-%d", ord);
        models[nm++] = mk(name, "poly", ord, 0, 0, th, p, c, w, kx);
    }
    int hasHistory = st->history_count >= 1;
    double xprev[2 * MMC_MAX_N];
    for (int i = 0; i < N; i++) {
        xprev[i] = st->prev2[i] * sc;
        xprev[N + i] = st->prev1[i] * sc;
    }
    if (hasHistory) {
        for (int eta = 0; eta <= 1; eta++) {
            for (int ord = 1; ord <= 2; ord++) {
                double c[2] = {0, 0}, w[2] = {2, 2}, th[2];
                fit_sample_pred(xn, xprev, 2 * N, N, ord, eta, c, w, th);
                char name[32];
                snprintf(name, sizeof(name), "samp.-%d-%d", ord, eta);
                models[nm++] = mk(name, "pred_samples", ord, eta, 0, th, ord, c, w, kx);
            }
        }
    }
    if (st->has_prev_model && st->prev_model.p > 0) {
        MmcModel *pm = &st->prev_model;
        if (strcmp(pm->family, "sin") == 0 || strcmp(pm->family, "poly") == 0 ||
            strcmp(pm->family, "pred_samples") == 0) {
            int factors[2] = {2, 10};
            for (int fi = 0; fi < 2; fi++) {
                int factor = factors[fi];
                double c[MMC_MAX_THETA], w[MMC_MAX_THETA], th[MMC_MAX_THETA];
                int p = pm->p;
                for (int j = 0; j < p; j++) {
                    c[j] = st->prev_theta_q[j];
                    w[j] = pm->theta_width[j] / factor;
                }
                if (strcmp(pm->family, "sin") == 0) fit_sinusoid(xn, N, cfg->fs, c, w, th);
                else if (strcmp(pm->family, "poly") == 0)
                    fit_poly(xn, N, pm->order, c, w, th);
                else
                    fit_sample_pred(xn, xprev, 2 * N, N, pm->order, pm->eta, c, w, th);
                char name[40];
                snprintf(name, sizeof(name), "para.-%d-%s", factor, pm->name);
                models[nm] = mk(name, "pred_para", pm->order, pm->eta, factor, th, p, c, w, kx);
                strncpy(models[nm].base_family, pm->family, sizeof(models[nm].base_family) - 1);
                nm++;
            }
        }
    }

    int M = nm;
    const char *rmethods[3] = {"none", "DCT", "DWT"};
    int Lr = 3;
    int nmb = (int)ceil(log2((double)M));
    if (nmb < 1) nmb = 1;
    int nl = (int)ceil(log2((double)Lr));
    if (nl < 1) nl = 1;

    int found = 0;
    double best_mse = 1e300;
    memset(best, 0, sizeof(*best));

    for (int im = 0; im < M; im++) {
        MmcModel *m = &models[im];
        int p = m->p;
        int nx_list[128];
        int nnx = 0;
        if (p == 0) nx_list[nnx++] = 0;
        else {
            int nmin = cfg->min_bits_theta * p;
            int nmax = cfg->max_bits_theta * p;
            if (nmax > cfg->n_tot) nmax = cfg->n_tot;
            for (int nx = nmin; nx <= nmax; nx += cfg->nx_step) nx_list[nnx++] = nx;
            if (nnx == 0) continue;
        }
        for (int ix = 0; ix < nnx; ix++) {
            int nx = nx_list[ix];
            int n_nx = 0;
            if (p > 0) {
                n_nx = (int)ceil(log2((double)(cfg->max_bits_theta * p + 1)));
                if (n_nx < 1) n_nx = 1;
            }
            int header_stage1 = nmb + cfg->n_kx + n_nx;
            if (header_stage1 + nx + nl > cfg->n_tot) continue;

            int bits[MMC_MAX_THETA];
            double theta_q[MMC_MAX_THETA];
            int theta_idx[MMC_MAX_THETA];
            double x_model[MMC_MAX_N];
            memset(bits, 0, sizeof(bits));
            memset(theta_q, 0, sizeof(theta_q));
            memset(theta_idx, 0, sizeof(theta_idx));
            if (p == 0) {
                memset(x_model, 0, (size_t)N * sizeof(double));
            } else {
                allocate_bits_greedy(m, nx, st, N, cfg->fs, cfg->min_bits_theta, cfg->max_bits_theta, bits);
                quantize_vector(m->theta_hat, bits, m->theta_center, m->theta_width, p, theta_q, theta_idx);
                reconstruct_model(m, theta_q, st, N, cfg->fs, x_model);
            }
            double r[MMC_MAX_N];
            for (int i = 0; i < N; i++) r[i] = xn[i] - x_model[i];

            for (int il = 0; il < Lr; il++) {
                const char *rname = rmethods[il];
                int nkr = 0, nnr = 0;
                if (strcmp(rname, "none") != 0) {
                    nkr = cfg->n_kr;
                    nnr = cfg->n_nr;
                }
                int nr_budget = cfg->n_tot - (header_stage1 + nx + nl + nkr + nnr);
                if (nr_budget < 0) continue;
                MmcResidual residual;
                encode_residual(r, N, rname, nr_budget, cfg->n_kr, 10, 1, &residual);
                double mse = 0.0;
                double xn_hat[MMC_MAX_N];
                for (int i = 0; i < N; i++) {
                    xn_hat[i] = x_model[i] + residual.rhat[i];
                    double d = xn[i] - xn_hat[i];
                    mse += d * d;
                }
                mse /= N;
                int total_bits = header_stage1 + nx + nl + nkr + nnr + residual.nbits;
                if (mse < best_mse) {
                    best_mse = mse;
                    found = 1;
                    memset(best, 0, sizeof(*best));
                    best->model = *m;
                    memcpy(best->theta_q, (p ? theta_q : theta_q), sizeof(best->theta_q));
                    best->residual = residual;
                    best->kx = kx;
                    best->nx = nx;
                    best->nr = residual.nbits;
                    best->header_bits = header_stage1 + nl + nkr + nnr;
                    best->total_bits = total_bits;
                    best->mse = mse;
                    best->nm = nmb;
                    best->nl = nl;
                    double p2 = pow(2.0, kx);
                    for (int i = 0; i < N; i++) best->xhat[i] = xn_hat[i] * p2;
                    memcpy(best->model.theta_hat, m->theta_hat, sizeof(m->theta_hat));
                }
            }
        }
    }
    free(xn);
    if (!found) return -1;
    memcpy(st->prev2, st->prev1, (size_t)N * sizeof(double));
    memcpy(st->prev1, best->xhat, (size_t)N * sizeof(double));
    st->prev_model = best->model;
    memcpy(st->prev_theta_q, best->theta_q, sizeof(st->prev_theta_q));
    st->has_prev_model = 1;
    st->history_count++;
    return 0;
}

int mmc_codec_run(const double *x, int n, const MmcOpts *opt, double *xhat, EwcrRunResult *out)
{
    memset(out, 0, sizeof(*out));
    strcpy(out->algo, "MMC");
    MmcOpts D;
    memset(&D, 0, sizeof(D));
    D.fs = 12800;
    D.fn = 50;
    D.N = 256;
    D.n_tot = 256;
    D.raw_bits_per_sample = 16;
    D.min_bits_theta = 2;
    D.max_bits_theta = 10;
    D.nx_step = 2;
    D.n_kx = 5;
    D.n_kr = 5;
    D.n_nr = 10;
    D.max_windows = 1000000;
    if (opt) D = *opt;
    int N = D.N;
    if (N > MMC_MAX_N) N = MMC_MAX_N;
    int W = n / N;
    if (W > D.max_windows) W = D.max_windows;
    if (W < 1) {
        snprintf(out->note, sizeof(out->note), "signal too short");
        return -1;
    }
    double t0 = ewcr_now_s();
    MmcState st;
    memset(&st, 0, sizeof(st));
    MmcPacket *pk = (MmcPacket *)ewcr_xcalloc((size_t)W, sizeof(MmcPacket));
    int bits_sum = 0;
    memset(xhat, 0, (size_t)n * sizeof(double));
    for (int w = 0; w < W; w++) {
        if (encode_window(x + w * N, N, &st, &D, &pk[w]) != 0) {
            snprintf(out->note, sizeof(out->note), "no feasible window %d", w);
            free(pk);
            return -1;
        }
        memcpy(xhat + w * N, pk[w].xhat, (size_t)N * sizeof(double));
        bits_sum += pk[w].total_bits;
    }
    double t1 = ewcr_now_s();

    /* independent decoder */
    MmcState dst;
    memset(&dst, 0, sizeof(dst));
    for (int w = 0; w < W; w++) {
        double xn_model[MMC_MAX_N], rn[MMC_MAX_N];
        reconstruct_model(&pk[w].model, pk[w].theta_q, &dst, N, D.fs, xn_model);
        decode_residual(&pk[w].residual, N, rn);
        double p2 = pow(2.0, pk[w].kx);
        for (int i = 0; i < N; i++) xhat[w * N + i] = (xn_model[i] + rn[i]) * p2;
        memcpy(dst.prev2, dst.prev1, (size_t)N * sizeof(double));
        memcpy(dst.prev1, xhat + w * N, (size_t)N * sizeof(double));
        dst.prev_model = pk[w].model;
        memcpy(dst.prev_theta_q, pk[w].theta_q, sizeof(dst.prev_theta_q));
        dst.has_prev_model = 1;
        dst.history_count++;
    }
    double t2 = ewcr_now_s();
    int nuse = W * N;
    ewcr_unified_metrics(x, xhat, nuse, (double)bits_sum, D.fs, D.fn, D.raw_bits_per_sample, &out->metrics);
    out->enc_s = t1 - t0;
    out->dec_s = t2 - t1;
    out->N = nuse;
    out->ok = 1;
    snprintf(out->note, sizeof(out->note), "windows=%d N=%d last=%s/%s", W, N, pk[W - 1].model.name,
             pk[W - 1].residual.name);
    free(pk);
    return 0;
}
