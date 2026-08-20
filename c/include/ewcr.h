#ifndef EWCR_H
#define EWCR_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double re, im;
} cpx;

static inline cpx cpx_make(double r, double i)
{
    cpx z;
    z.re = r;
    z.im = i;
    return z;
}
static inline cpx cpx_add(cpx a, cpx b) { return cpx_make(a.re + b.re, a.im + b.im); }
static inline cpx cpx_sub(cpx a, cpx b) { return cpx_make(a.re - b.re, a.im - b.im); }
static inline cpx cpx_mul(cpx a, cpx b)
{
    return cpx_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline cpx cpx_scale(cpx a, double s) { return cpx_make(a.re * s, a.im * s); }
static inline cpx cpx_conj(cpx a) { return cpx_make(a.re, -a.im); }
static inline double cpx_abs2(cpx a) { return a.re * a.re + a.im * a.im; }
static inline double cpx_abs(cpx a) { return sqrt(cpx_abs2(a)); }

/* ---- memory / time / rng ------------------------------------------------ */
void *ewcr_xmalloc(size_t n);
void *ewcr_xcalloc(size_t n, size_t sz);
void ewcr_free(void *p);
double ewcr_now_s(void);

typedef struct {
    uint64_t s;
} EwcrRng;
void ewcr_rng_seed(EwcrRng *r, uint64_t seed);
double ewcr_rng_u01(EwcrRng *r);
double ewcr_rng_randn(EwcrRng *r);
int ewcr_rng_randi(EwcrRng *r, int lo, int hi);

/* ---- math kernels ------------------------------------------------------- */
int ewcr_nextpow2(int n);
void ewcr_fft(cpx *x, int n, int inverse);          /* MATLAB fft/ifft convention */
void ewcr_fft_real(const double *x, int n, int fft_len, cpx *X);
void ewcr_interpft(const cpx *x, int n, cpx *y, int ny);
void ewcr_dct_matrix(double *T, int n);              /* orthonormal DCT-II, row-major */
void ewcr_dct_fwd(const double *x, double *c, int n);
void ewcr_dct_inv(const double *c, double *x, int n);
void ewcr_haar_fwd(const double *x, double *c, int n);
void ewcr_haar_inv(const double *c, double *x, int n);
int ewcr_ls_solve(const double *A, int m, int n, const double *b, double *x);
int ewcr_ls_solve_cpx(const cpx *A, int m, int n, const cpx *b, cpx *x);

void ewcr_db4_filters(double *lod, double *hid, double *lor, double *hir);
int ewcr_wavedec_db4(const double *x, int n, int L, double *C, int *book, int *ncoef);
int ewcr_waverec_db4(const double *C, const int *book, int L, double *x);

int ewcr_huff_encode(const int32_t *sig, int n, uint8_t **bits, int *nbits,
                     int32_t **syms, int *nsyms, uint8_t **codebits, int **codelen);
int ewcr_huff_decode(const uint8_t *bits, int nbits, const int32_t *syms, int nsyms,
                     const uint8_t *codebits, const int *codelen, int32_t *out, int nout);

/* ---- metrics / signals -------------------------------------------------- */
typedef struct {
    int N;
    double bits_original;
    double bits_compressed;
    double CR;
    double rate_bps;
    double NMSE_lin;
    double NMSE_dB;
    double SNR_dB;
    double RMSE;
    double PRD;
    double MAXE;
    double PSNR_dB;
    double fund_amp;
    double fund_amp_hat;
    double fund_amp_relerr;
    double fund_phase_err_deg;
    double THD;
    double THD_hat;
    double THD_abs_error;
} EwcrMetrics;

void ewcr_unified_metrics(const double *x, const double *xhat, int n,
                          double bits_compressed, double fs, double f0,
                          int raw_bits, EwcrMetrics *m);

int ewcr_gen_ieee1159(const char *kind, double fs, int cycles, double *x, int *n_out);

typedef struct {
    char path[512];
    double fs;
    int n;
    double *UA, *IA, *UB, *IB, *UC, *IC;
} EwcrWaveCsv;

int ewcr_load_wave_csv(const char *path, EwcrWaveCsv *w);
void ewcr_free_wave_csv(EwcrWaveCsv *w);
const double *ewcr_wave_channel(const EwcrWaveCsv *w, const char *name);

/* ---- codec options ------------------------------------------------------ */
typedef struct {
    double W, act_thresh, energy_frac;
    int Rk, block_cycles, ke_max, RQ;
} AsbcOpts;

typedef struct {
    int level;          /* <=0 => auto */
    char thr_mode[16];  /* "universal" or "energy" */
    double thr_frac;
    int qbits;
    int delta;
    int RQ;
} DwtOpts;

typedef struct {
    double M_ratio, tol;
    int K0, seed, ybits, RQ;
} CsOpts;

typedef struct {
    double fs, fn;
    int N, n_tot, raw_bits_per_sample;
    int min_bits_theta, max_bits_theta, nx_step;
    int n_kx, n_kr, n_nr;
    int verbose;
    int max_windows;
} MmcOpts;

typedef struct {
    double fs, f0, signal_scale, G, beta;
    double input_full_scale_pu, fft_full_scale_pu;
    int Nppc, lagrange_order;
    int input_quantize, input_bits, fft_quantize, fft_bits;
    int original_bits_per_sample, use_lzw, verbose;
    char sync_mode[16];
} SvdcsOpts;

typedef struct {
    double fs, f0;
    int raw_bits_per_sample;
    AsbcOpts asbc;
    DwtOpts dwt;
    CsOpts cs;
    MmcOpts mmc;
    SvdcsOpts svdcs;
} EwcrBenchCfg;

void ewcr_default_cfg(EwcrBenchCfg *cfg);

typedef struct {
    char algo[32];
    char signal[32];
    int N;
    int ok;
    double enc_s, dec_s;
    EwcrMetrics metrics;
    char note[160];
} EwcrRunResult;

int asbc_codec_run(const double *x, int n, double fs, const AsbcOpts *opt,
                   double *xhat, EwcrRunResult *out);
int dwt_codec_run(const double *x, int n, const DwtOpts *opt,
                  double *xhat, EwcrRunResult *out);
int cs_omp_codec_run(const double *x, int n, const CsOpts *opt,
                     double *xhat, EwcrRunResult *out);
int mmc_codec_run(const double *x, int n, const MmcOpts *opt,
                  double *xhat, EwcrRunResult *out);
int svdcs_codec_run(const double *x, int n, const SvdcsOpts *opt,
                    double *xhat, EwcrRunResult *out);

int ewcr_run_codec(const char *algo, const double *x, int n, const EwcrBenchCfg *cfg,
                   double *xhat, EwcrRunResult *out);

int ewcr_unit_tests(int verbose);

#ifdef __cplusplus
}
#endif

#endif /* EWCR_H */
