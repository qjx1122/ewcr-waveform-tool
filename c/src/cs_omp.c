#include "ewcr.h"

#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double vec_norm(const cpx *v, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += cpx_abs2(v[i]);
    return sqrt(s);
}

static void omp_recover(const cpx *A, int M, int N, const cpx *y, int K0, double tol, cpx *s_hat)
{
    memset(s_hat, 0, (size_t)N * sizeof(cpx));
    cpx *r = (cpx *)ewcr_xmalloc((size_t)M * sizeof(cpx));
    memcpy(r, y, (size_t)M * sizeof(cpx));
    int *supp = (int *)ewcr_xcalloc((size_t)K0, sizeof(int));
    cpx *As = (cpx *)ewcr_xcalloc((size_t)M * (size_t)K0, sizeof(cpx));
    cpx *coef = (cpx *)ewcr_xcalloc((size_t)K0, sizeof(cpx));
    double target = tol * vec_norm(y, M);
    int k = 0;
    while (k < K0 && vec_norm(r, M) > target) {
        int idx = 0;
        double best = -1.0;
        for (int j = 0; j < N; j++) {
            cpx acc = cpx_make(0, 0);
            for (int i = 0; i < M; i++) acc = cpx_add(acc, cpx_mul(cpx_conj(A[i * N + j]), r[i]));
            double a = cpx_abs(acc);
            if (a > best) {
                best = a;
                idx = j;
            }
        }
        for (int t = 0; t < k; t++) {
            if (idx == supp[t]) {
                /* zero that correlation by excluding: pick next best not in support */
                best = -1.0;
                idx = 0;
                for (int j = 0; j < N; j++) {
                    int used = 0;
                    for (int u = 0; u < k; u++) if (supp[u] == j) used = 1;
                    if (used) continue;
                    cpx acc = cpx_make(0, 0);
                    for (int i = 0; i < M; i++) acc = cpx_add(acc, cpx_mul(cpx_conj(A[i * N + j]), r[i]));
                    double a = cpx_abs(acc);
                    if (a > best) {
                        best = a;
                        idx = j;
                    }
                }
                break;
            }
        }
        supp[k] = idx;
        for (int i = 0; i < M; i++) As[i * K0 + k] = A[i * N + idx];
        k++;
        /* least squares on As[:,0:k] * coef = y. As stored with lda=K0 */
        cpx *As_pack = (cpx *)ewcr_xmalloc((size_t)M * (size_t)k * sizeof(cpx));
        for (int i = 0; i < M; i++)
            for (int j = 0; j < k; j++) As_pack[i * k + j] = As[i * K0 + j];
        ewcr_ls_solve_cpx(As_pack, M, k, y, coef);
        for (int i = 0; i < M; i++) {
            cpx acc = cpx_make(0, 0);
            for (int j = 0; j < k; j++) acc = cpx_add(acc, cpx_mul(As[i * K0 + j], coef[j]));
            r[i] = cpx_sub(y[i], acc);
        }
        free(As_pack);
    }
    if (k > 0) {
        cpx *As_pack = (cpx *)ewcr_xmalloc((size_t)M * (size_t)k * sizeof(cpx));
        for (int i = 0; i < M; i++)
            for (int j = 0; j < k; j++) As_pack[i * k + j] = As[i * K0 + j];
        ewcr_ls_solve_cpx(As_pack, M, k, y, coef);
        for (int j = 0; j < k; j++) s_hat[supp[j]] = coef[j];
        free(As_pack);
    }
    free(r);
    free(supp);
    free(As);
    free(coef);
}

static void form_A_dft(const double *Phi, int M, int N, cpx *A)
{
    /* A[m,k] = DFT_k(Phi[m,:]) / sqrt(N) */
    double s = 1.0 / sqrt((double)N);
    cpx *row = (cpx *)ewcr_xmalloc((size_t)N * sizeof(cpx));
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) row[n] = cpx_make(Phi[m * N + n], 0);
        ewcr_fft(row, N, 0);
        for (int k = 0; k < N; k++) A[m * N + k] = cpx_scale(row[k], s);
    }
    free(row);
}

static void psi_times_s(const cpx *s, int N, double *xhat)
{
    cpx *tmp = (cpx *)ewcr_xmalloc((size_t)N * sizeof(cpx));
    memcpy(tmp, s, (size_t)N * sizeof(cpx));
    ewcr_fft(tmp, N, 0);
    double sc = 1.0 / sqrt((double)N);
    for (int i = 0; i < N; i++) xhat[i] = tmp[i].re * sc;
    free(tmp);
}

int cs_omp_codec_run(const double *x, int n, const CsOpts *opt, double *xhat, EwcrRunResult *out)
{
    memset(out, 0, sizeof(*out));
    strcpy(out->algo, "CS-OMP");
    CsOpts D;
    memset(&D, 0, sizeof(D));
    D.M_ratio = 0.20;
    D.tol = 1e-6;
    D.K0 = 48;
    D.seed = 42;
    D.ybits = 16;
    D.RQ = 16;
    if (opt) D = *opt;
    double t0 = ewcr_now_s();

    int N = n;
    int M = (int)round(D.M_ratio * N);
    if (M < 1) M = 1;
    if (M > N) M = N;

    EwcrRng rng;
    ewcr_rng_seed(&rng, (uint64_t)D.seed);
    double *Phi = (double *)ewcr_xmalloc((size_t)M * (size_t)N * sizeof(double));
    double ism = 1.0 / sqrt((double)M);
    for (int i = 0; i < M * N; i++) Phi[i] = ewcr_rng_randn(&rng) * ism;

    cpx *A = (cpx *)ewcr_xmalloc((size_t)M * (size_t)N * sizeof(cpx));
    form_A_dft(Phi, M, N, A);

    cpx *y = (cpx *)ewcr_xcalloc((size_t)M, sizeof(cpx));
    for (int m = 0; m < M; m++) {
        double s = 0.0;
        const double *row = Phi + (size_t)m * N;
        for (int i = 0; i < N; i++) s += row[i] * x[i];
        y[m].re = s;
    }

    double yscale = 0.0;
    for (int m = 0; m < M; m++) {
        double a = fabs(y[m].re);
        if (a > yscale) yscale = a;
    }
    if (yscale == 0.0) yscale = 1.0;
    int32_t *yq = (int32_t *)ewcr_xmalloc((size_t)M * sizeof(int32_t));
    double qmax = (double)(1 << (D.ybits - 1));
    for (int m = 0; m < M; m++) {
        double v = round(y[m].re / yscale * qmax);
        if (v < -qmax) v = -qmax;
        if (v > qmax - 1) v = qmax - 1;
        yq[m] = (int32_t)v;
    }
    double bits_total = (double)M * D.ybits + 160.0;
    double t1 = ewcr_now_s();

    /* decode: dequant + OMP */
    cpx *ydec = (cpx *)ewcr_xcalloc((size_t)M, sizeof(cpx));
    for (int m = 0; m < M; m++) ydec[m].re = (double)yq[m] * yscale / qmax;

    cpx *s_hat = (cpx *)ewcr_xcalloc((size_t)N, sizeof(cpx));
    omp_recover(A, M, N, ydec, D.K0, D.tol, s_hat);
    psi_times_s(s_hat, N, xhat);
    double t2 = ewcr_now_s();

    double fs = (D.fs > 0.0) ? D.fs : 12800.0;
    double f0 = (D.f0 > 0.0) ? D.f0 : 50.0;
    ewcr_unified_metrics(x, xhat, n, bits_total, fs, f0, D.RQ, &out->metrics);
    out->enc_s = t1 - t0;
    out->dec_s = t2 - t1;
    out->N = n;
    out->ok = 1;
    snprintf(out->note, sizeof(out->note), "M=%d N=%d K0=%d", M, N, D.K0);

    {
        EwcrBuf b;
        ewcr_buf_init(&b);
        ewcr_pack_header(&b, 3, N, fs, f0, bits_total);
        ewcr_buf_u32(&b, (uint32_t)M);
        ewcr_buf_u32(&b, (uint32_t)D.K0);
        ewcr_buf_i32(&b, D.seed);
        ewcr_buf_u32(&b, (uint32_t)D.ybits);
        ewcr_buf_f64(&b, D.tol);
        ewcr_buf_f64(&b, yscale);
        for (int i = 0; i < M; i++) ewcr_buf_i32(&b, yq[i]);
        out->compressed = b.d;
        out->compressed_nbytes = (int)b.n;
    }

    free(Phi);
    free(A);
    free(y);
    free(yq);
    free(ydec);
    free(s_hat);
    return 0;
}
