#include "ewcr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void make_seg(int N, double start_frac, double dur_frac, double *seg)
{
    int a = (int)round(start_frac * N);
    if (a < 1) a = 1;
    int b = a + (int)round(dur_frac * N);
    if (b > N) b = N;
    memset(seg, 0, (size_t)N * sizeof(double));
    for (int i = a - 1; i < b; i++) seg[i] = 1.0; /* 1-based MATLAB a:b */
}

int ewcr_gen_ieee1159(const char *kind, double fs, int cycles, double *x, int *n_out)
{
    const double F0 = 50.0;
    int N = (int)round(fs * cycles / F0);
    if (N < 8) return -1;
    *n_out = N;
    double w = 2.0 * M_PI * F0 / fs;

    double sag_mag = 0.7, sag_start = 0.2, sag_dur = 0.3;
    double swell_mag = 1.3, swell_start = 0.2, swell_dur = 0.3;
    double harm_a[4] = {0.05, 0.03, 0.02, 0.01};
    int harm_h[4] = {3, 5, 7, 11};
    double trans_amp = 0.6, trans_f = 550, trans_start = 0.35;
    double notch_depth = 0.25, notch_frac = 0.05;
    double flicker_m = 0.05, flicker_f = 8;
    double spike_amp = 1.2;
    int spike_n = 3;

    double *base = (double *)ewcr_xmalloc((size_t)N * sizeof(double));
    double *seg = (double *)ewcr_xmalloc((size_t)N * sizeof(double));
    for (int n = 0; n < N; n++) base[n] = sin(w * n);

    if (strcmp(kind, "pure") == 0) {
        memcpy(x, base, (size_t)N * sizeof(double));
    } else if (strcmp(kind, "sag") == 0) {
        make_seg(N, sag_start, sag_dur, seg);
        for (int i = 0; i < N; i++) x[i] = base[i] * (1.0 - (1.0 - sag_mag) * seg[i]);
    } else if (strcmp(kind, "swell") == 0) {
        make_seg(N, swell_start, swell_dur, seg);
        for (int i = 0; i < N; i++) x[i] = base[i] * (1.0 + (swell_mag - 1.0) * seg[i]);
    } else if (strcmp(kind, "interruption") == 0) {
        make_seg(N, 0.2, 0.3, seg);
        for (int i = 0; i < N; i++) x[i] = base[i] * (1.0 - seg[i]);
    } else if (strcmp(kind, "harmonics") == 0) {
        memcpy(x, base, (size_t)N * sizeof(double));
        double sa = 0.0;
        for (int i = 0; i < 4; i++) {
            sa += harm_a[i];
            for (int n = 0; n < N; n++) x[n] += harm_a[i] * sin(harm_h[i] * w * n + 0.3 * (i + 1));
        }
        for (int n = 0; n < N; n++) x[n] /= (1.0 + sa);
    } else if (strcmp(kind, "osc_transient") == 0) {
        memcpy(x, base, (size_t)N * sizeof(double));
        make_seg(N, trans_start, 0.05, seg);
        for (int n = 0; n < N; n++) {
            double env = exp(-(n - trans_start * N) / (0.02 * N));
            x[n] += trans_amp * sin(2 * M_PI * trans_f / fs * n) * env * seg[n];
        }
    } else if (strcmp(kind, "notch") == 0) {
        memcpy(x, base, (size_t)N * sizeof(double));
        int period = (int)round(fs / F0);
        int notch_w = (int)round(notch_frac * period);
        if (notch_w < 1) notch_w = 1;
        for (double cyc = 0.25; cyc <= cycles - 0.75 + 1e-9; cyc += 1.0) {
            int p = (int)round(period * cyc); /* 1-based peak_pos */
            if (p > 0 && p <= N) {
                int lo = p - notch_w;
                if (lo < 1) lo = 1;
                int hi = p + notch_w;
                if (hi > N) hi = N;
                for (int i = lo - 1; i < hi; i++) x[i] = (1.0 - notch_depth) * x[i];
            }
        }
    } else if (strcmp(kind, "flicker") == 0) {
        for (int n = 0; n < N; n++)
            x[n] = base[n] * (1.0 + flicker_m * sin(2 * M_PI * flicker_f / fs * n));
    } else if (strcmp(kind, "spike") == 0) {
        memcpy(x, base, (size_t)N * sizeof(double));
        EwcrRng rng;
        ewcr_rng_seed(&rng, 7);
        for (int i = 0; i < spike_n; i++) {
            double t = (spike_n == 1) ? 0.5 : (0.1 + 0.8 * i / (double)(spike_n - 1));
            int p = (int)round(t * N) + ewcr_rng_randi(&rng, 1, 5);
            if (p < 2) p = 2;
            if (p > N - 1) p = N - 1;
            x[p - 1] = spike_amp * ((x[p - 1] >= 0) ? 1.0 : -1.0);
        }
    } else if (strcmp(kind, "complex") == 0) {
        memcpy(x, base, (size_t)N * sizeof(double));
        make_seg(N, 0.15, 0.4, seg);
        for (int i = 0; i < N; i++) x[i] *= (1.0 - (1.0 - 0.75) * seg[i]);
        double sa = 0.0;
        for (int i = 0; i < 4; i++) {
            sa += harm_a[i];
            for (int n = 0; n < N; n++) x[n] += harm_a[i] * sin(harm_h[i] * w * n + 0.3 * (i + 1));
        }
        for (int n = 0; n < N; n++) x[n] /= (1.0 + sa);
        make_seg(N, 0.55, 0.02, seg);
        for (int n = 0; n < N; n++) {
            double env = exp(-(n - 0.55 * N) / (0.01 * N));
            x[n] += 0.5 * sin(2 * M_PI * 700 / fs * n) * env * seg[n];
        }
    } else {
        free(base);
        free(seg);
        return -2;
    }
    free(base);
    free(seg);
    return 0;
}

static double parse_tod_seconds(const char *ts)
{
    const char *p = strrchr(ts, ' ');
    p = p ? p + 1 : ts;
    int h = 0, m = 0;
    double s = 0.0;
    if (sscanf(p, "%d:%d:%lf", &h, &m, &s) != 3) return NAN;
    return (double)h * 3600.0 + (double)m * 60.0 + s;
}

int ewcr_load_wave_csv(const char *path, EwcrWaveCsv *w)
{
    memset(w, 0, sizeof(*w));
    strncpy(w->path, path, sizeof(w->path) - 1);
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    char line[1024];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return -2;
    }
    int cap = 4096;
    int n = 0;
    double *ch[6];
    double *tod = (double *)ewcr_xmalloc((size_t)cap * sizeof(double));
    for (int c = 0; c < 6; c++) ch[c] = (double *)ewcr_xmalloc((size_t)cap * sizeof(double));
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '\0' || line[0] == '\n' || line[0] == '\r') continue;
        char *tok = strtok(line, ",");
        if (!tok) continue;
        double t = parse_tod_seconds(tok);
        double v[6];
        int got = 0;
        for (int c = 0; c < 6; c++) {
            tok = strtok(NULL, ",");
            if (!tok || tok[0] == '\0' || tok[0] == '\n') {
                v[c] = 0.0;
            } else {
                v[c] = strtod(tok, NULL);
                got++;
            }
        }
        if (got < 1) continue;
        if (n >= cap) {
            cap *= 2;
            tod = (double *)realloc(tod, (size_t)cap * sizeof(double));
            for (int c = 0; c < 6; c++) ch[c] = (double *)realloc(ch[c], (size_t)cap * sizeof(double));
            if (!tod) {
                fprintf(stderr, "csv realloc failed\n");
                exit(2);
            }
        }
        tod[n] = t;
        for (int c = 0; c < 6; c++) ch[c][n] = v[c];
        n++;
    }
    fclose(fp);
    if (n < 8) {
        free(tod);
        for (int c = 0; c < 6; c++) free(ch[c]);
        return -3;
    }
    double acc = 0.0;
    int ndt = 0;
    for (int i = 1; i < n && i < 2000; i++) {
        double d = tod[i] - tod[i - 1];
        if (d > 1e-9 && d < 0.1) {
            acc += d;
            ndt++;
        }
    }
    free(tod);
    double dt = (ndt > 0) ? (acc / ndt) : 1e-4;
    w->fs = 1.0 / dt;
    w->n = n;
    w->UA = ch[0];
    w->IA = ch[1];
    w->UB = ch[2];
    w->IB = ch[3];
    w->UC = ch[4];
    w->IC = ch[5];
    return 0;
}

void ewcr_free_wave_csv(EwcrWaveCsv *w)
{
    if (!w) return;
    free(w->UA);
    free(w->IA);
    free(w->UB);
    free(w->IB);
    free(w->UC);
    free(w->IC);
    memset(w, 0, sizeof(*w));
}

const double *ewcr_wave_channel(const EwcrWaveCsv *w, const char *name)
{
    if (!w || !name) return NULL;
    if (strcmp(name, "UA") == 0) return w->UA;
    if (strcmp(name, "IA") == 0) return w->IA;
    if (strcmp(name, "UB") == 0) return w->UB;
    if (strcmp(name, "IB") == 0) return w->IB;
    if (strcmp(name, "UC") == 0) return w->UC;
    if (strcmp(name, "IC") == 0) return w->IC;
    return NULL;
}
