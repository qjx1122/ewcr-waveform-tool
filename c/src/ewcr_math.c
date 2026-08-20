#include "ewcr.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void *ewcr_xmalloc(size_t n)
{
    if (n == 0) n = 1;
    void *p = malloc(n);
    if (!p) {
        fprintf(stderr, "ewcr: out of memory (%zu bytes)\n", n);
        exit(2);
    }
    return p;
}

void *ewcr_xcalloc(size_t n, size_t sz)
{
    if (n == 0) n = 1;
    if (sz == 0) sz = 1;
    void *p = calloc(n, sz);
    if (!p) {
        fprintf(stderr, "ewcr: out of memory (%zu x %zu)\n", n, sz);
        exit(2);
    }
    return p;
}

void ewcr_free(void *p) { free(p); }

double ewcr_now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* SplitMix64 -> xorshift for portable deterministic RNG (not MATLAB mt19937). */
void ewcr_rng_seed(EwcrRng *r, uint64_t seed)
{
    seed += 0x9E3779B97F4A7C15ULL;
    seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9ULL;
    seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EBULL;
    r->s = seed ^ (seed >> 31);
    if (r->s == 0) r->s = 0xA5A5A5A5A5A5A5A5ULL;
}

double ewcr_rng_u01(EwcrRng *r)
{
    uint64_t x = r->s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    r->s = x;
    uint64_t u = x * 0x2545F4914F6CDD1DULL;
    return (double)(u >> 11) * (1.0 / 9007199254740992.0);
}

double ewcr_rng_randn(EwcrRng *r)
{
    double u, v, s;
    do {
        u = 2.0 * ewcr_rng_u01(r) - 1.0;
        v = 2.0 * ewcr_rng_u01(r) - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    return u * sqrt(-2.0 * log(s) / s);
}

int ewcr_rng_randi(EwcrRng *r, int lo, int hi)
{
    if (hi < lo) {
        int t = lo;
        lo = hi;
        hi = t;
    }
    int span = hi - lo + 1;
    return lo + (int)floor(ewcr_rng_u01(r) * span);
}

int ewcr_nextpow2(int n)
{
    if (n <= 1) return 0;
    int p = 0;
    int v = 1;
    while (v < n) {
        v <<= 1;
        p++;
        if (p > 30) break;
    }
    return p;
}

static int is_pow2(int n) { return n > 0 && (n & (n - 1)) == 0; }

static void bit_reverse(cpx *x, int n)
{
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            cpx t = x[i];
            x[i] = x[j];
            x[j] = t;
        }
    }
}

static void fft_pow2(cpx *x, int n, int inverse)
{
    bit_reverse(x, n);
    for (int len = 2; len <= n; len <<= 1) {
        double ang = (inverse ? 2.0 : -2.0) * M_PI / (double)len;
        cpx wlen = cpx_make(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cpx w = cpx_make(1.0, 0.0);
            int half = len >> 1;
            for (int j = 0; j < half; j++) {
                cpx u = x[i + j];
                cpx v = cpx_mul(x[i + j + half], w);
                x[i + j] = cpx_add(u, v);
                x[i + j + half] = cpx_sub(u, v);
                w = cpx_mul(w, wlen);
            }
        }
    }
    if (inverse) {
        double invn = 1.0 / (double)n;
        for (int i = 0; i < n; i++) x[i] = cpx_scale(x[i], invn);
    }
}

static void fft_bluestein(cpx *x, int n, int inverse)
{
    int m = 1 << ewcr_nextpow2(2 * n - 1);
    cpx *a = (cpx *)ewcr_xcalloc((size_t)m, sizeof(cpx));
    cpx *b = (cpx *)ewcr_xcalloc((size_t)m, sizeof(cpx));
    cpx *w = (cpx *)ewcr_xmalloc((size_t)n * sizeof(cpx));
    double sign = inverse ? 1.0 : -1.0;
    for (int k = 0; k < n; k++) {
        double ang = sign * M_PI * (double)k * (double)k / (double)n;
        w[k] = cpx_make(cos(ang), sin(ang));
        a[k] = cpx_mul(x[k], w[k]);
        b[k] = cpx_make(cos(-ang), sin(-ang));
    }
    for (int k = 1; k < n; k++) b[m - k] = b[k];
    fft_pow2(a, m, 0);
    fft_pow2(b, m, 0);
    for (int k = 0; k < m; k++) a[k] = cpx_mul(a[k], b[k]);
    fft_pow2(a, m, 1);
    for (int k = 0; k < n; k++) x[k] = cpx_mul(a[k], w[k]);
    if (inverse) {
        double invn = 1.0 / (double)n;
        for (int k = 0; k < n; k++) x[k] = cpx_scale(x[k], invn);
    }
    free(a);
    free(b);
    free(w);
}

void ewcr_fft(cpx *x, int n, int inverse)
{
    if (n <= 1) return;
    if (is_pow2(n)) fft_pow2(x, n, inverse);
    else fft_bluestein(x, n, inverse);
}

void ewcr_fft_real(const double *x, int n, int fft_len, cpx *X)
{
    for (int i = 0; i < fft_len; i++) {
        X[i].re = (i < n) ? x[i] : 0.0;
        X[i].im = 0.0;
    }
    ewcr_fft(X, fft_len, 0);
}

void ewcr_interpft(const cpx *x, int n, cpx *y, int ny)
{
    if (n <= 0 || ny <= 0) return;
    if (n == 1) {
        for (int i = 0; i < ny; i++) y[i] = x[0];
        return;
    }
    if (n == ny) {
        memcpy(y, x, (size_t)n * sizeof(cpx));
        return;
    }
    cpx *X = (cpx *)ewcr_xmalloc((size_t)n * sizeof(cpx));
    memcpy(X, x, (size_t)n * sizeof(cpx));
    ewcr_fft(X, n, 0);
    cpx *Y = (cpx *)ewcr_xcalloc((size_t)ny, sizeof(cpx));
    int nyqst = (n + 1) / 2; /* ceil((n+1)/2) with 0-based copy count = nyqst */
    /* MATLAB: nyqst = ceil((m+1)/2) is 1-based length of positive incl DC */
    for (int i = 0; i < nyqst; i++) Y[i] = X[i];
    int tail = n - nyqst;
    for (int i = 0; i < tail; i++) Y[ny - tail + i] = X[nyqst + i];
    if ((n % 2) == 0) {
        /* split Nyquist bin */
        cpx half = cpx_scale(Y[nyqst - 1], 0.5);
        Y[nyqst - 1] = half;
        Y[ny - n + nyqst - 1] = half;
    }
    ewcr_fft(Y, ny, 1);
    double scale = (double)ny / (double)n;
    /* ifft already divided by ny; MATLAB does ifft(Y)*(ny/m).
       Our ifft divides by ny, MATLAB ifft divides by ny then * ny/m = 1/m * sum.
       MATLAB: y = ifft(Y)* (ny/m). ifft already /ny, so extra * ny/m => overall /m.
       Our ifft /ny, we need overall /n = /m, so multiply by ny/n. Yes. */
    for (int i = 0; i < ny; i++) y[i] = cpx_scale(Y[i], scale);
    free(X);
    free(Y);
}

void ewcr_dct_matrix(double *T, int n)
{
    for (int k = 0; k < n; k++) {
        double a = (k == 0) ? sqrt(1.0 / n) : sqrt(2.0 / n);
        for (int i = 0; i < n; i++) {
            T[k * n + i] = a * cos(M_PI * (0.5 + i) * k / n);
        }
    }
}

void ewcr_dct_fwd(const double *x, double *c, int n)
{
    double *T = (double *)ewcr_xmalloc((size_t)n * (size_t)n * sizeof(double));
    ewcr_dct_matrix(T, n);
    for (int k = 0; k < n; k++) {
        double s = 0.0;
        const double *row = T + (size_t)k * n;
        for (int i = 0; i < n; i++) s += row[i] * x[i];
        c[k] = s;
    }
    free(T);
}

void ewcr_dct_inv(const double *c, double *x, int n)
{
    double *T = (double *)ewcr_xmalloc((size_t)n * (size_t)n * sizeof(double));
    ewcr_dct_matrix(T, n);
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int k = 0; k < n; k++) s += T[k * n + i] * c[k];
        x[i] = s;
    }
    free(T);
}

void ewcr_haar_fwd(const double *x, double *c, int n)
{
    double *temp = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
    memcpy(temp, x, (size_t)n * sizeof(double));
    memset(c, 0, (size_t)n * sizeof(double));
    int len = n;
    int pos = n;
    const double s2 = 1.0 / sqrt(2.0);
    while (len > 1) {
        int n2 = len / 2;
        double *a = (double *)ewcr_xmalloc((size_t)n2 * sizeof(double));
        double *d = (double *)ewcr_xmalloc((size_t)n2 * sizeof(double));
        for (int i = 0; i < n2; i++) {
            a[i] = (temp[2 * i] + temp[2 * i + 1]) * s2;
            d[i] = (temp[2 * i] - temp[2 * i + 1]) * s2;
        }
        memcpy(c + (pos - n2), d, (size_t)n2 * sizeof(double));
        pos -= n2;
        memcpy(temp, a, (size_t)n2 * sizeof(double));
        len = n2;
        free(a);
        free(d);
    }
    c[0] = temp[0];
    free(temp);
}

void ewcr_haar_inv(const double *c, double *x, int n)
{
    const double s2 = 1.0 / sqrt(2.0);
    double *a = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
    a[0] = c[0];
    int pos = 1;
    int len = 1;
    while (len < n) {
        double *temp = (double *)ewcr_xcalloc((size_t)(2 * len), sizeof(double));
        for (int i = 0; i < len; i++) {
            double d = c[pos + i];
            temp[2 * i] = (a[i] + d) * s2;
            temp[2 * i + 1] = (a[i] - d) * s2;
        }
        pos += len;
        memcpy(a, temp, (size_t)(2 * len) * sizeof(double));
        free(temp);
        len *= 2;
    }
    memcpy(x, a, (size_t)n * sizeof(double));
    free(a);
}

/* Column-pivoted normal equations, A is m x n row-major. */
int ewcr_ls_solve(const double *A, int m, int n, const double *b, double *x)
{
    if (n <= 0) return 0;
    double *G = (double *)ewcr_xcalloc((size_t)n * (size_t)n, sizeof(double));
    double *g = (double *)ewcr_xcalloc((size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double s = 0.0;
            for (int k = 0; k < m; k++) s += A[k * n + i] * A[k * n + j];
            G[i * n + j] = G[j * n + i] = s;
        }
        double t = 0.0;
        for (int k = 0; k < m; k++) t += A[k * n + i] * b[k];
        g[i] = t;
    }
    double tr = 0.0;
    for (int i = 0; i < n; i++) tr += G[i * n + i];
    double ridge = 1e-14 * (tr > 1.0 ? tr : 1.0);
    for (int i = 0; i < n; i++) G[i * n + i] += ridge;

    /* Gaussian elimination with partial pivot */
    int *piv = (int *)ewcr_xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) piv[i] = i;
    double *M = G;
    for (int k = 0; k < n; k++) {
        int best = k;
        double bv = fabs(M[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(M[i * n + k]);
            if (v > bv) {
                bv = v;
                best = i;
            }
        }
        if (bv < 1e-18) {
            memset(x, 0, (size_t)n * sizeof(double));
            free(G);
            free(g);
            free(piv);
            return -1;
        }
        if (best != k) {
            for (int j = 0; j < n; j++) {
                double tmp = M[k * n + j];
                M[k * n + j] = M[best * n + j];
                M[best * n + j] = tmp;
            }
            double tg = g[k];
            g[k] = g[best];
            g[best] = tg;
        }
        double akk = M[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double f = M[i * n + k] / akk;
            for (int j = k; j < n; j++) M[i * n + j] -= f * M[k * n + j];
            g[i] -= f * g[k];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = g[i];
        for (int j = i + 1; j < n; j++) s -= M[i * n + j] * x[j];
        x[i] = s / M[i * n + i];
    }
    free(G);
    free(g);
    free(piv);
    return 0;
}

int ewcr_ls_solve_cpx(const cpx *A, int m, int n, const cpx *b, cpx *x)
{
    /* Solve real 2n system from normal equations A^H A x = A^H b */
    if (n <= 0) return 0;
    int N = n;
    double *G = (double *)ewcr_xcalloc((size_t)N * N, sizeof(double)); /* real part of Hermitian */
    double *Gi = (double *)ewcr_xcalloc((size_t)N * N, sizeof(double));
    cpx *rhs = (cpx *)ewcr_xcalloc((size_t)N, sizeof(cpx));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cpx s = cpx_make(0, 0);
            for (int k = 0; k < m; k++) {
                s = cpx_add(s, cpx_mul(cpx_conj(A[k * N + i]), A[k * N + j]));
            }
            G[i * N + j] = s.re;
            Gi[i * N + j] = s.im;
        }
        cpx t = cpx_make(0, 0);
        for (int k = 0; k < m; k++) t = cpx_add(t, cpx_mul(cpx_conj(A[k * N + i]), b[k]));
        rhs[i] = t;
    }
    /* Build 2N real system [G, -Gi; Gi, G] [xr; xi] = [re; im] and solve it
     * directly (G is already A^H A). Do NOT form normal equations again. */
    int n2 = 2 * N;
    double *R = (double *)ewcr_xcalloc((size_t)n2 * n2, sizeof(double));
    double *rb = (double *)ewcr_xcalloc((size_t)n2, sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * n2 + j] = G[i * N + j];
            R[i * n2 + (j + N)] = -Gi[i * N + j];
            R[(i + N) * n2 + j] = Gi[i * N + j];
            R[(i + N) * n2 + (j + N)] = G[i * N + j];
        }
        rb[i] = rhs[i].re;
        rb[i + N] = rhs[i].im;
    }
    double tr = 0.0;
    for (int i = 0; i < n2; i++) tr += R[i * n2 + i];
    double ridge = 1e-14 * (tr > 1.0 ? tr : 1.0);
    for (int i = 0; i < n2; i++) R[i * n2 + i] += ridge;
    for (int k = 0; k < n2; k++) {
        int best = k;
        double bv = fabs(R[k * n2 + k]);
        for (int i = k + 1; i < n2; i++) {
            double v = fabs(R[i * n2 + k]);
            if (v > bv) {
                bv = v;
                best = i;
            }
        }
        if (bv < 1e-18) {
            memset(x, 0, (size_t)N * sizeof(cpx));
            free(G);
            free(Gi);
            free(rhs);
            free(R);
            free(rb);
            return -1;
        }
        if (best != k) {
            for (int j = 0; j < n2; j++) {
                double tmp = R[k * n2 + j];
                R[k * n2 + j] = R[best * n2 + j];
                R[best * n2 + j] = tmp;
            }
            double tg = rb[k];
            rb[k] = rb[best];
            rb[best] = tg;
        }
        double akk = R[k * n2 + k];
        for (int i = k + 1; i < n2; i++) {
            double f = R[i * n2 + k] / akk;
            for (int j = k; j < n2; j++) R[i * n2 + j] -= f * R[k * n2 + j];
            rb[i] -= f * rb[k];
        }
    }
    double *rx = (double *)ewcr_xcalloc((size_t)n2, sizeof(double));
    for (int i = n2 - 1; i >= 0; i--) {
        double s = rb[i];
        for (int j = i + 1; j < n2; j++) s -= R[i * n2 + j] * rx[j];
        rx[i] = s / R[i * n2 + i];
    }
    for (int i = 0; i < N; i++) x[i] = cpx_make(rx[i], rx[i + N]);
    free(G);
    free(Gi);
    free(rhs);
    free(R);
    free(rb);
    free(rx);
    return 0;
}

void ewcr_db4_filters(double *lod, double *hid, double *lor, double *hir)
{
    /* Orthonormal db4 taps (same analysis/synthesis pair, periodized PR). */
    static const double h[8] = {
        0.23037781330885523,
        0.7148465705525415,
        0.6308807679298587,
        -0.027983769416983849,
        -0.18703481171888114,
        0.030841381835987667,
        0.032883011666982945,
        -0.010597401784997278
    };
    const int L = 8;
    for (int k = 0; k < L; k++) {
        lod[k] = h[k];
        lor[k] = h[k];
        hid[k] = ((k & 1) ? -1.0 : 1.0) * h[L - 1 - k];
        hir[k] = hid[k];
    }
}

static void dwt1_per(const double *x, int n, const double *lod, const double *hid, int L,
                     double *a, double *d)
{
    int n2 = n / 2;
    for (int i = 0; i < n2; i++) {
        double sa = 0.0, sd = 0.0;
        for (int k = 0; k < L; k++) {
            int idx = (2 * i - k) % n;
            if (idx < 0) idx += n;
            sa += lod[k] * x[idx];
            sd += hid[k] * x[idx];
        }
        a[i] = sa;
        d[i] = sd;
    }
}

static void idwt1_per(const double *a, const double *d, int n2, const double *lor, const double *hir,
                      int L, double *x)
{
    int n = n2 * 2;
    memset(x, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < n2; i++) {
        for (int k = 0; k < L; k++) {
            int idx = (2 * i - k) % n;
            if (idx < 0) idx += n;
            x[idx] += lor[k] * a[i] + hir[k] * d[i];
        }
    }
}

int ewcr_wavedec_db4(const double *x, int n, int L, double *C, int *book, int *ncoef)
{
    double lod[8], hid[8], lor[8], hir[8];
    ewcr_db4_filters(lod, hid, lor, hir);
    if (L < 1) L = 1;
    int need = n;
    for (int i = 0; i < L; i++) {
        if (need % 2) return -1;
        need /= 2;
    }
    double *cur = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
    memcpy(cur, x, (size_t)n * sizeof(double));
    int cur_n = n;
    double **details = (double **)ewcr_xmalloc((size_t)L * sizeof(double *));
    int *dlen = (int *)ewcr_xmalloc((size_t)L * sizeof(int));
    for (int lev = 0; lev < L; lev++) {
        int n2 = cur_n / 2;
        double *a = (double *)ewcr_xmalloc((size_t)n2 * sizeof(double));
        double *d = (double *)ewcr_xmalloc((size_t)n2 * sizeof(double));
        dwt1_per(cur, cur_n, lod, hid, 8, a, d);
        details[lev] = d; /* level 1 stored first; MATLAB book is cD_L ... cD_1 */
        dlen[lev] = n2;
        free(cur);
        cur = a;
        cur_n = n2;
    }
    /* Pack MATLAB-style: cA_L, cD_L, ..., cD_1 */
    int pos = 0;
    book[0] = cur_n;
    memcpy(C + pos, cur, (size_t)cur_n * sizeof(double));
    pos += cur_n;
    for (int lev = L - 1; lev >= 0; lev--) {
        book[L - lev] = dlen[lev];
        memcpy(C + pos, details[lev], (size_t)dlen[lev] * sizeof(double));
        pos += dlen[lev];
        free(details[lev]);
    }
    book[L + 1] = n;
    *ncoef = pos;
    free(cur);
    free(details);
    free(dlen);
    (void)lor;
    (void)hir;
    return 0;
}

int ewcr_waverec_db4(const double *C, const int *book, int L, double *x)
{
    double lod[8], hid[8], lor[8], hir[8];
    ewcr_db4_filters(lod, hid, lor, hir);
    int nA = book[0];
    double *a = (double *)ewcr_xmalloc((size_t)book[L + 1] * sizeof(double));
    memcpy(a, C, (size_t)nA * sizeof(double));
    int a_n = nA;
    int pos = nA;
    for (int j = 0; j < L; j++) {
        int nD = book[j + 1];
        const double *d = C + pos;
        double *rec = (double *)ewcr_xmalloc((size_t)(a_n * 2) * sizeof(double));
        idwt1_per(a, d, a_n, lor, hir, 8, rec);
        free(a);
        a = rec;
        a_n *= 2;
        pos += nD;
        (void)nD;
    }
    memcpy(x, a, (size_t)a_n * sizeof(double));
    free(a);
    (void)lod;
    (void)hid;
    return 0;
}

/* ---------------- Huffman ------------------------------------------------ */
typedef struct HuffNode {
    int32_t sym;
    int freq;
    int is_leaf;
    struct HuffNode *lo, *hi;
} HuffNode;

static HuffNode *huff_new(int32_t sym, int freq, int is_leaf, HuffNode *lo, HuffNode *hi)
{
    HuffNode *n = (HuffNode *)ewcr_xcalloc(1, sizeof(HuffNode));
    n->sym = sym;
    n->freq = freq;
    n->is_leaf = is_leaf;
    n->lo = lo;
    n->hi = hi;
    return n;
}

static void huff_free(HuffNode *n)
{
    if (!n) return;
    huff_free(n->lo);
    huff_free(n->hi);
    free(n);
}

static void huff_walk(HuffNode *n, uint8_t *buf, int depth,
                      int32_t *syms, uint8_t *codebits, int *codelen, int *ns)
{
    if (!n) return;
    if (n->is_leaf) {
        int i = *ns;
        syms[i] = n->sym;
        codelen[i] = depth > 0 ? depth : 1;
        uint8_t acc = 0;
        /* store first 8 bits only for tiny codes; full codes kept in codebits as packed length<=32 */
        uint32_t full = 0;
        for (int b = 0; b < depth; b++) full = (full << 1) | (buf[b] & 1);
        if (depth == 0) {
            codebits[i] = 0;
            codelen[i] = 1;
        } else {
            /* we store only via parallel arrays later; keep lowest 8 in codebits placeholder */
            (void)acc;
            memcpy(&codebits[i], &full, 1); /* placeholder, overwritten by int table */
        }
        (void)full;
        (*ns)++;
        return;
    }
    buf[depth] = 0;
    huff_walk(n->lo, buf, depth + 1, syms, codebits, codelen, ns);
    buf[depth] = 1;
    huff_walk(n->hi, buf, depth + 1, syms, codebits, codelen, ns);
}

typedef struct {
    int32_t sym;
    uint32_t code;
    int len;
} HuffCode;

static void huff_fill(HuffNode *n, uint32_t code, int depth, HuffCode *tab, int *ns)
{
    if (!n) return;
    if (n->is_leaf) {
        tab[*ns].sym = n->sym;
        tab[*ns].code = (depth == 0) ? 0 : code;
        tab[*ns].len = depth == 0 ? 1 : depth;
        (*ns)++;
        return;
    }
    huff_fill(n->lo, code << 1, depth + 1, tab, ns);
    huff_fill(n->hi, (code << 1) | 1, depth + 1, tab, ns);
}

int ewcr_huff_encode(const int32_t *sig, int n, uint8_t **bits, int *nbits,
                     int32_t **syms_out, int *nsyms, uint8_t **codebits, int **codelen)
{
    *bits = NULL;
    *nbits = 0;
    *syms_out = NULL;
    *nsyms = 0;
    *codebits = NULL;
    *codelen = NULL;
    if (n <= 0) return 0;

    /* unique symbols */
    int32_t *tmp = (int32_t *)ewcr_xmalloc((size_t)n * sizeof(int32_t));
    memcpy(tmp, sig, (size_t)n * sizeof(int32_t));
    /* insertion sort unique */
    for (int i = 1; i < n; i++) {
        int32_t v = tmp[i];
        int j = i - 1;
        while (j >= 0 && tmp[j] > v) {
            tmp[j + 1] = tmp[j];
            j--;
        }
        tmp[j + 1] = v;
    }
    int nu = 1;
    for (int i = 1; i < n; i++) if (tmp[i] != tmp[nu - 1]) tmp[nu++] = tmp[i];
    int *freq = (int *)ewcr_xcalloc((size_t)nu, sizeof(int));
    for (int i = 0; i < n; i++) {
        int lo = 0, hi = nu - 1, mid = 0;
        while (lo <= hi) {
            mid = (lo + hi) / 2;
            if (tmp[mid] < sig[i]) lo = mid + 1;
            else if (tmp[mid] > sig[i]) hi = mid - 1;
            else break;
        }
        freq[mid]++;
    }

    HuffCode *tab = (HuffCode *)ewcr_xcalloc((size_t)nu, sizeof(HuffCode));
    int ntab = 0;
    if (nu == 1) {
        tab[0].sym = tmp[0];
        tab[0].code = 0;
        tab[0].len = 1;
        ntab = 1;
    } else {
        HuffNode **heap = (HuffNode **)ewcr_xmalloc((size_t)nu * sizeof(HuffNode *));
        int hs = nu;
        for (int i = 0; i < nu; i++) heap[i] = huff_new(tmp[i], freq[i], 1, NULL, NULL);
        while (hs > 1) {
            int i0 = 0, i1 = 1;
            if (heap[i1]->freq < heap[i0]->freq) {
                int t = i0;
                i0 = i1;
                i1 = t;
            }
            for (int i = 2; i < hs; i++) {
                if (heap[i]->freq < heap[i0]->freq) {
                    i1 = i0;
                    i0 = i;
                } else if (heap[i]->freq < heap[i1]->freq) i1 = i;
            }
            HuffNode *a = heap[i0], *b = heap[i1];
            HuffNode *p = huff_new(0, a->freq + b->freq, 0, a, b);
            int ia = i0 > i1 ? i0 : i1;
            int ib = i0 < i1 ? i0 : i1;
            heap[ia] = heap[hs - 1];
            hs--;
            heap[ib] = heap[hs - 1];
            hs--;
            heap[hs++] = p;
        }
        huff_fill(heap[0], 0, 0, tab, &ntab);
        huff_free(heap[0]);
        free(heap);
    }

    *nsyms = ntab;
    *syms_out = (int32_t *)ewcr_xmalloc((size_t)ntab * sizeof(int32_t));
    *codelen = (int *)ewcr_xmalloc((size_t)ntab * sizeof(int));
    *codebits = (uint8_t *)ewcr_xcalloc((size_t)ntab * 4, 1); /* unused; codes in parallel int */
    /* stash codes in extra int array via codebits pointer? keep local map */
    uint32_t *codes = (uint32_t *)ewcr_xmalloc((size_t)ntab * sizeof(uint32_t));
    for (int i = 0; i < ntab; i++) {
        (*syms_out)[i] = tab[i].sym;
        (*codelen)[i] = tab[i].len;
        codes[i] = tab[i].code;
        /* pack 32-bit code into 4 bytes of codebits */
    }
    free(*codebits);
    *codebits = (uint8_t *)codes;

    int cap = n * 16 + 64;
    uint8_t *bitbuf = (uint8_t *)ewcr_xcalloc((size_t)cap, 1);
    int nb = 0;
    for (int i = 0; i < n; i++) {
        int k = -1;
        for (int j = 0; j < ntab; j++) if ((*syms_out)[j] == sig[i]) {
            k = j;
            break;
        }
        if (k < 0) continue;
        uint32_t c = codes[k];
        int L = (*codelen)[k];
        for (int b = L - 1; b >= 0; b--) {
            if (nb >= cap) {
                cap *= 2;
                bitbuf = (uint8_t *)realloc(bitbuf, (size_t)cap);
                if (!bitbuf) {
                    fprintf(stderr, "huff realloc failed\n");
                    exit(2);
                }
            }
            bitbuf[nb++] = (uint8_t)((c >> b) & 1);
        }
    }
    *bits = bitbuf;
    *nbits = nb;
    free(tmp);
    free(freq);
    free(tab);
    return 0;
}

int ewcr_huff_decode(const uint8_t *bits, int nbits, const int32_t *syms, int nsyms,
                     const uint8_t *codebits, const int *codelen, int32_t *out, int nout)
{
    if (nout <= 0) return 0;
    if (nsyms <= 0) {
        memset(out, 0, (size_t)nout * sizeof(int32_t));
        return 0;
    }
    if (nsyms == 1) {
        for (int i = 0; i < nout; i++) out[i] = syms[0];
        return 0;
    }
    const uint32_t *codes = (const uint32_t *)codebits;
    int pos = 0;
    int produced = 0;
    while (produced < nout && pos < nbits) {
        uint32_t acc = 0;
        int matched = 0;
        for (int len = 1; len <= 32 && pos < nbits; len++) {
            acc = (acc << 1) | (uint32_t)(bits[pos++] & 1);
            for (int s = 0; s < nsyms; s++) {
                if (codelen[s] == len && codes[s] == acc) {
                    out[produced++] = syms[s];
                    matched = 1;
                    break;
                }
            }
            if (matched) break;
        }
        if (!matched) break;
    }
    while (produced < nout) out[produced++] = 0;
    return 0;
}

void ewcr_unified_metrics(const double *x, const double *xhat, int n,
                          double bits_compressed, double fs, double f0,
                          int raw_bits, EwcrMetrics *m)
{
    memset(m, 0, sizeof(*m));
    m->N = n;
    m->bits_original = (double)n * raw_bits;
    m->bits_compressed = bits_compressed;
    m->CR = m->bits_original / fmax(bits_compressed, 1.0);
    m->rate_bps = bits_compressed / fmax((double)n, 1.0);
    double Es = 0.0, Ee = 0.0, peak = 0.0;
    for (int i = 0; i < n; i++) {
        double e = x[i] - xhat[i];
        Es += x[i] * x[i];
        Ee += e * e;
        double ax = fabs(x[i]);
        if (ax > peak) peak = ax;
        double ae = fabs(e);
        if (ae > m->MAXE) m->MAXE = ae;
    }
    m->NMSE_lin = Ee / fmax(Es, 1e-30);
    m->NMSE_dB = 10.0 * log10(fmax(m->NMSE_lin, 1e-300));
    m->SNR_dB = -m->NMSE_dB;
    m->RMSE = sqrt(Ee / fmax((double)n, 1.0));
    m->PRD = 100.0 * sqrt(Ee / fmax(Es, 1e-30));
    m->PSNR_dB = 20.0 * log10(fmax(peak, 1e-30) / fmax(m->RMSE, 1e-30));

    /* Pearson r(x, xhat): waveform similarity vs reconstruction. */
    {
        double mx = 0.0, mh = 0.0;
        for (int i = 0; i < n; i++) {
            mx += x[i];
            mh += xhat[i];
        }
        mx /= fmax((double)n, 1.0);
        mh /= fmax((double)n, 1.0);
        double num = 0.0, vx = 0.0, vh = 0.0;
        for (int i = 0; i < n; i++) {
            double dx = x[i] - mx;
            double dh = xhat[i] - mh;
            num += dx * dh;
            vx += dx * dx;
            vh += dh * dh;
        }
        double den = sqrt(vx * vh);
        m->corr = (den > 1e-30) ? (num / den) : 0.0;
        if (m->corr > 1.0) m->corr = 1.0;
        if (m->corr < -1.0) m->corr = -1.0;
        m->similarity_pct = 100.0 * m->corr;
    }

    if (fs <= 0 || f0 <= 0) return;
    double *A = (double *)ewcr_xmalloc((size_t)n * 2 * sizeof(double));
    double *co = (double *)ewcr_xmalloc(2 * sizeof(double));
    double *ch = (double *)ewcr_xmalloc(2 * sizeof(double));
    for (int i = 0; i < n; i++) {
        double t = (double)i / fs;
        A[i * 2 + 0] = cos(2 * M_PI * f0 * t);
        A[i * 2 + 1] = sin(2 * M_PI * f0 * t);
    }
    ewcr_ls_solve(A, n, 2, x, co);
    ewcr_ls_solve(A, n, 2, xhat, ch);
    m->fund_amp = hypot(co[0], co[1]);
    m->fund_amp_hat = hypot(ch[0], ch[1]);
    m->fund_amp_relerr = fabs(m->fund_amp_hat - m->fund_amp) / fmax(m->fund_amp, 1e-30);
    double pho = atan2(-co[1], co[0]);
    double phh = atan2(-ch[1], ch[0]);
    double dp = phh - pho;
    dp = atan2(sin(dp), cos(dp));
    m->fund_phase_err_deg = fabs(dp) * 180.0 / M_PI;

    int Hmax = (int)floor((fs / 2.0) / f0);
    if (Hmax > 25) Hmax = 25;
    if (Hmax < 1) Hmax = 1;
    double *Ah = (double *)ewcr_xcalloc((size_t)Hmax, sizeof(double));
    double *Ahh = (double *)ewcr_xcalloc((size_t)Hmax, sizeof(double));
    for (int h = 1; h <= Hmax; h++) {
        for (int i = 0; i < n; i++) {
            double t = (double)i / fs;
            A[i * 2 + 0] = cos(2 * M_PI * h * f0 * t);
            A[i * 2 + 1] = sin(2 * M_PI * h * f0 * t);
        }
        ewcr_ls_solve(A, n, 2, x, co);
        ewcr_ls_solve(A, n, 2, xhat, ch);
        Ah[h - 1] = hypot(co[0], co[1]);
        Ahh[h - 1] = hypot(ch[0], ch[1]);
    }
    double sh = 0.0, shh = 0.0;
    for (int h = 1; h < Hmax; h++) {
        sh += Ah[h] * Ah[h];
        shh += Ahh[h] * Ahh[h];
    }
    m->THD = sqrt(sh) / fmax(Ah[0], 1e-30);
    m->THD_hat = sqrt(shh) / fmax(Ahh[0], 1e-30);
    m->THD_abs_error = fabs(m->THD_hat - m->THD);
    free(A);
    free(co);
    free(ch);
    free(Ah);
    free(Ahh);
}

void ewcr_default_cfg(EwcrBenchCfg *cfg)
{
    memset(cfg, 0, sizeof(*cfg));
    cfg->fs = 12800;
    cfg->f0 = 50;
    cfg->raw_bits_per_sample = 16;

    cfg->asbc.W = 20;
    cfg->asbc.Rk = 8;
    cfg->asbc.block_cycles = 5;
    cfg->asbc.act_thresh = 1e-4;
    cfg->asbc.ke_max = 48;
    cfg->asbc.energy_frac = 0.99;
    cfg->asbc.RQ = 16;

    cfg->dwt.level = 0;
    strcpy(cfg->dwt.thr_mode, "energy");
    cfg->dwt.thr_frac = 0.99;
    cfg->dwt.qbits = 8;
    cfg->dwt.delta = 1;
    cfg->dwt.RQ = 16;

    cfg->cs.M_ratio = 0.20;
    cfg->cs.K0 = 48;
    cfg->cs.tol = 1e-6;
    cfg->cs.seed = 42;
    cfg->cs.ybits = 16;
    cfg->cs.RQ = 16;

    cfg->mmc.fs = cfg->fs;
    cfg->mmc.fn = cfg->f0;
    cfg->mmc.N = (int)round(cfg->fs / cfg->f0);
    cfg->mmc.n_tot = cfg->mmc.N;
    cfg->mmc.raw_bits_per_sample = 16;
    cfg->mmc.min_bits_theta = 2;
    cfg->mmc.max_bits_theta = 10;
    cfg->mmc.nx_step = 2;
    cfg->mmc.n_kx = 5;
    cfg->mmc.n_kr = 5;
    cfg->mmc.n_nr = 10;
    cfg->mmc.verbose = 0;
    cfg->mmc.max_windows = 1000000;

    cfg->svdcs.fs = cfg->fs;
    cfg->svdcs.f0 = cfg->f0;
    cfg->svdcs.Nppc = (int)round(cfg->fs / cfg->f0);
    cfg->svdcs.signal_scale = 1.0;
    cfg->svdcs.G = 0.03;
    cfg->svdcs.beta = 0.10;
    cfg->svdcs.input_full_scale_pu = 2.0;
    cfg->svdcs.fft_full_scale_pu = 2.0;
    cfg->svdcs.lagrange_order = 3;
    cfg->svdcs.input_quantize = 1;
    cfg->svdcs.input_bits = 16;
    cfg->svdcs.fft_quantize = 1;
    cfg->svdcs.fft_bits = 16;
    cfg->svdcs.original_bits_per_sample = 16;
    cfg->svdcs.use_lzw = 1;
    cfg->svdcs.verbose = 0;
    strcpy(cfg->svdcs.sync_mode, "nominal");
}

int ewcr_run_codec(const char *algo, const double *x, int n, const EwcrBenchCfg *cfg,
                   double *xhat, EwcrRunResult *out)
{
    memset(out, 0, sizeof(*out));
    strncpy(out->algo, algo, sizeof(out->algo) - 1);
    if (strcmp(algo, "ASBC") == 0)
        return asbc_codec_run(x, n, cfg->fs, &cfg->asbc, xhat, out);
    if (strcmp(algo, "DWT-Hybrid") == 0 || strcmp(algo, "DWT") == 0)
        return dwt_codec_run(x, n, &cfg->dwt, xhat, out);
    if (strcmp(algo, "CS-OMP") == 0 || strcmp(algo, "CS") == 0)
        return cs_omp_codec_run(x, n, &cfg->cs, xhat, out);
    if (strcmp(algo, "MMC") == 0)
        return mmc_codec_run(x, n, &cfg->mmc, xhat, out);
    if (strcmp(algo, "SVDCS") == 0)
        return svdcs_codec_run(x, n, &cfg->svdcs, xhat, out);
    snprintf(out->note, sizeof(out->note), "unknown codec");
    return -1;
}

int ewcr_unit_tests(int verbose)
{
    int fails = 0;
#define CHECK(cond, msg)                                                                 \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            fprintf(stderr, "UNIT FAIL: %s\n", msg);                                     \
            fails++;                                                                     \
        } else if (verbose) {                                                            \
            printf("  PASS %s\n", msg);                                                  \
        }                                                                                \
    } while (0)

    /* FFT round-trip */
    {
        int n = 128;
        cpx *x = (cpx *)ewcr_xmalloc((size_t)n * sizeof(cpx));
        cpx *y = (cpx *)ewcr_xmalloc((size_t)n * sizeof(cpx));
        for (int i = 0; i < n; i++) x[i] = cpx_make(sin(2 * M_PI * 3 * i / (double)n), 0.1 * i);
        memcpy(y, x, (size_t)n * sizeof(cpx));
        ewcr_fft(y, n, 0);
        ewcr_fft(y, n, 1);
        double e = 0;
        for (int i = 0; i < n; i++) e += cpx_abs2(cpx_sub(x[i], y[i]));
        CHECK(sqrt(e / n) < 1e-12, "fft pow2 roundtrip");
        free(x);
        free(y);
    }
    {
        int n = 100;
        cpx *x = (cpx *)ewcr_xmalloc((size_t)n * sizeof(cpx));
        cpx *y = (cpx *)ewcr_xmalloc((size_t)n * sizeof(cpx));
        for (int i = 0; i < n; i++) x[i] = cpx_make(cos(0.3 * i), sin(0.2 * i));
        memcpy(y, x, (size_t)n * sizeof(cpx));
        ewcr_fft(y, n, 0);
        ewcr_fft(y, n, 1);
        double e = 0;
        for (int i = 0; i < n; i++) e += cpx_abs2(cpx_sub(x[i], y[i]));
        CHECK(sqrt(e / n) < 1e-10, "fft bluestein roundtrip");
        free(x);
        free(y);
    }
    /* Haar / DCT PR */
    {
        int n = 128;
        double *x = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
        double *c = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
        double *y = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
        for (int i = 0; i < n; i++) x[i] = sin(2 * M_PI * i / n) + 0.2 * cos(6 * M_PI * i / n);
        ewcr_haar_fwd(x, c, n);
        ewcr_haar_inv(c, y, n);
        double e = 0;
        for (int i = 0; i < n; i++) e += (x[i] - y[i]) * (x[i] - y[i]);
        CHECK(sqrt(e / n) < 1e-12, "haar roundtrip");
        ewcr_dct_fwd(x, c, n);
        ewcr_dct_inv(c, y, n);
        e = 0;
        for (int i = 0; i < n; i++) e += (x[i] - y[i]) * (x[i] - y[i]);
        CHECK(sqrt(e / n) < 1e-11, "dct roundtrip");
        free(x);
        free(c);
        free(y);
    }
    /* db4 PR */
    {
        int n = 256, L = 4, ncoef = 0;
        double *x = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
        double *C = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
        double *y = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
        int book[16];
        for (int i = 0; i < n; i++) x[i] = sin(2 * M_PI * 5 * i / (double)n);
        int rc = ewcr_wavedec_db4(x, n, L, C, book, &ncoef);
        CHECK(rc == 0, "wavedec rc");
        ewcr_waverec_db4(C, book, L, y);
        double e = 0;
        for (int i = 0; i < n; i++) e += (x[i] - y[i]) * (x[i] - y[i]);
        CHECK(sqrt(e / n) < 1e-11, "db4 wavedec/waverec PR");
        free(x);
        free(C);
        free(y);
    }
    /* Huffman */
    {
        int32_t sig[16] = {1, 1, 2, 3, 3, 3, 0, -1, -1, 2, 2, 2, 7, 7, 7, 7};
        uint8_t *bits = NULL, *cb = NULL;
        int32_t *syms = NULL;
        int *cl = NULL, nbits = 0, nsyms = 0;
        ewcr_huff_encode(sig, 16, &bits, &nbits, &syms, &nsyms, &cb, &cl);
        int32_t out[16];
        ewcr_huff_decode(bits, nbits, syms, nsyms, cb, cl, out, 16);
        int ok = 1;
        for (int i = 0; i < 16; i++) if (out[i] != sig[i]) ok = 0;
        CHECK(ok, "huffman roundtrip");
        free(bits);
        free(cb);
        free(syms);
        free(cl);
    }
    /* LS */
    {
        double A[6] = {1, 0, 0, 1, 1, 1};
        double b[3] = {1, 2, 3};
        double x[2];
        ewcr_ls_solve(A, 3, 2, b, x);
        CHECK(fabs(x[0] - 1) < 1e-9 && fabs(x[1] - 2) < 1e-9, "ls 2-col");
    }
    /* Pearson identity / scaled copy */
    {
        double x[8], y[8];
        for (int i = 0; i < 8; i++) {
            x[i] = (double)i;
            y[i] = 2.0 * x[i] + 3.0;
        }
        EwcrMetrics m;
        ewcr_unified_metrics(x, x, 8, 8, 0, 0, 16, &m);
        CHECK(fabs(m.corr - 1.0) < 1e-12 && fabs(m.similarity_pct - 100.0) < 1e-9, "corr identity");
        ewcr_unified_metrics(x, y, 8, 8, 0, 0, 16, &m);
        CHECK(fabs(m.corr - 1.0) < 1e-12, "corr affine");
    }
#undef CHECK
    return fails;
}
