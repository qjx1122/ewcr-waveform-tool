#include "ewcr.h"

#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void quantize_signed(const double *x, int n, int nbits, double fullScale, int32_t *q, double *step)
{
    double maxCode = (double)((1 << (nbits - 1)) - 1);
    double minCode = -(double)(1 << (nbits - 1));
    *step = fullScale / maxCode;
    for (int i = 0; i < n; i++) {
        double v = round(x[i] / (*step));
        if (v > maxCode) v = maxCode;
        if (v < minCode) v = minCode;
        q[i] = (int32_t)v;
    }
}

static double lagrange_one(const double *y, int N, double t, int order)
{
    int m = order + 1;
    int left = (int)floor(t) - (int)floor((m - 1) / 2.0);
    /* MATLAB 1-based floor(t) */
    int t1 = (int)floor(t); /* 1-based sample index floor */
    left = t1 - (int)floor((m - 1) / 2.0);
    if (left < 1) left = 1;
    if (left > N - m + 1) left = N - m + 1;
    double nodes[8], vals[8], L[8];
    if (m > 8) m = 8;
    for (int a = 0; a < m; a++) {
        nodes[a] = (double)(left + a);
        vals[a] = y[left + a - 1];
        L[a] = 1.0;
    }
    double dmin = 1e300;
    int jj = 0;
    for (int a = 0; a < m; a++) {
        double d = fabs(nodes[a] - t);
        if (d < dmin) {
            dmin = d;
            jj = a;
        }
    }
    if (dmin < 10.0 * 2.22e-16 * fmax(1.0, fabs(t))) return vals[jj];
    for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) if (a != b)
            L[a] *= (t - nodes[b]) / (nodes[a] - nodes[b]);
    }
    double s = 0.0;
    for (int a = 0; a < m; a++) s += vals[a] * L[a];
    return s;
}

static int prepare_cycles_nominal(const double *x, int n, const SvdcsOpts *cfg, double *frames, int *nframes)
{
    double T = cfg->fs / cfg->f0;
    int M = cfg->Nppc;
    int nF = (int)floor(n / T);
    while (nF > 0) {
        double lastStart = 1.0 + (nF - 1) * T;
        double lastQuery = lastStart + (M - 1) * (T / M);
        if (lastQuery <= n + 10.0 * 2.22e-16 * n) break;
        nF--;
    }
    if (nF < 1) return -1;
    for (int k = 0; k < nF; k++) {
        double start = 1.0 + k * T;
        for (int i = 0; i < M; i++) {
            double tq = start + i * (T / M);
            frames[k * M + i] = lagrange_one(x, n, tq, cfg->lagrange_order);
        }
    }
    *nframes = nF;
    return 0;
}

/* 12-bit LZW */
static uint8_t *pack12(const uint16_t *codes, int n, int *nbytes)
{
    int nb = (int)ceil(12.0 * n / 8.0);
    uint8_t *packed = (uint8_t *)ewcr_xcalloc((size_t)(nb + 2), 1);
    int p = 0, i = 0;
    while (i + 1 < n) {
        uint16_t a = codes[i], b = codes[i + 1];
        packed[p] = (uint8_t)(a >> 4);
        packed[p + 1] = (uint8_t)(((a & 15) << 4) | (b >> 8));
        packed[p + 2] = (uint8_t)(b & 255);
        p += 3;
        i += 2;
    }
    if (i < n) {
        uint16_t a = codes[i];
        packed[p] = (uint8_t)(a >> 4);
        packed[p + 1] = (uint8_t)((a & 15) << 4);
    }
    *nbytes = nb;
    return packed;
}

static void unpack12(const uint8_t *packed, int nCodes, uint16_t *codes)
{
    int p = 0, i = 0;
    while (i + 1 < nCodes) {
        uint16_t b1 = packed[p], b2 = packed[p + 1], b3 = packed[p + 2];
        codes[i] = (uint16_t)((b1 << 4) | (b2 >> 4));
        codes[i + 1] = (uint16_t)(((b2 & 15) << 8) | b3);
        p += 3;
        i += 2;
    }
    if (i < nCodes) {
        uint16_t b1 = packed[p], b2 = packed[p + 1];
        codes[i] = (uint16_t)((b1 << 4) | (b2 >> 4));
    }
}

static uint8_t *lzw_encode12(const uint8_t *bytes, int nbytes, int *npacked, int *nCodes)
{
    if (nbytes <= 0) {
        *npacked = 0;
        *nCodes = 0;
        return (uint8_t *)ewcr_xmalloc(1);
    }
    const int MAXC = 4095;
    uint16_t *transition = (uint16_t *)ewcr_xcalloc((size_t)(MAXC + 1) * 256, sizeof(uint16_t));
    uint16_t nextCode = 256;
    uint16_t *codes = (uint16_t *)ewcr_xmalloc((size_t)(nbytes + 8) * sizeof(uint16_t));
    int nc = 0;
    uint16_t w = bytes[0];
    for (int i = 1; i < nbytes; i++) {
        uint16_t b = bytes[i];
        uint16_t z = transition[(int)w * 256 + b];
        if (z != 0) w = z;
        else {
            codes[nc++] = w;
            if (nextCode <= MAXC) {
                transition[(int)w * 256 + b] = nextCode;
                nextCode++;
            }
            w = b;
        }
    }
    codes[nc++] = w;
    uint8_t *packed = pack12(codes, nc, npacked);
    *nCodes = nc;
    free(transition);
    free(codes);
    return packed;
}

static uint8_t *lzw_decode12(const uint8_t *packed, int nCodes, int *nout)
{
    if (nCodes == 0) {
        *nout = 0;
        return (uint8_t *)ewcr_xmalloc(1);
    }
    uint16_t *codes = (uint16_t *)ewcr_xmalloc((size_t)nCodes * sizeof(uint16_t));
    unpack12(packed, nCodes, codes);
    uint8_t **dict = (uint8_t **)ewcr_xcalloc(4096, sizeof(uint8_t *));
    int *dlen = (int *)ewcr_xcalloc(4096, sizeof(int));
    for (int k = 0; k < 256; k++) {
        dict[k] = (uint8_t *)ewcr_xmalloc(1);
        dict[k][0] = (uint8_t)k;
        dlen[k] = 1;
    }
    int nextCode = 256;
    int cap = nCodes * 8 + 16;
    uint8_t *out = (uint8_t *)ewcr_xmalloc((size_t)cap);
    int noutb = 0;
    int oldCode = codes[0];
    memcpy(out, dict[oldCode], (size_t)dlen[oldCode]);
    noutb = dlen[oldCode];
    uint8_t *prev = (uint8_t *)ewcr_xmalloc((size_t)(cap));
    int prevn = dlen[oldCode];
    memcpy(prev, dict[oldCode], (size_t)prevn);
    for (int i = 1; i < nCodes; i++) {
        int c = codes[i];
        uint8_t *entry;
        int elen;
        uint8_t tmp[8];
        if (c < nextCode && dict[c]) {
            entry = dict[c];
            elen = dlen[c];
        } else if (c == nextCode) {
            if (noutb + prevn + 1 > cap) {
                cap *= 2;
                out = (uint8_t *)realloc(out, (size_t)cap);
                prev = (uint8_t *)realloc(prev, (size_t)cap);
            }
            memcpy(tmp, prev, 1); /* placeholder not used */
            /* entry = prev + prev[0] — allocate temp */
            uint8_t *e2 = (uint8_t *)ewcr_xmalloc((size_t)prevn + 1);
            memcpy(e2, prev, (size_t)prevn);
            e2[prevn] = prev[0];
            entry = e2;
            elen = prevn + 1;
            /* take ownership below */
            if (noutb + elen > cap) {
                cap = noutb + elen + 1024;
                out = (uint8_t *)realloc(out, (size_t)cap);
            }
            memcpy(out + noutb, entry, (size_t)elen);
            noutb += elen;
            if (nextCode <= 4095) {
                dict[nextCode] = (uint8_t *)ewcr_xmalloc((size_t)prevn + 1);
                memcpy(dict[nextCode], prev, (size_t)prevn);
                dict[nextCode][prevn] = entry[0];
                dlen[nextCode] = prevn + 1;
                nextCode++;
            }
            memcpy(prev, entry, (size_t)elen);
            prevn = elen;
            free(e2);
            continue;
        } else {
            fprintf(stderr, "Invalid LZW code %d\n", c);
            break;
        }
        if (noutb + elen > cap) {
            cap = (noutb + elen) * 2;
            out = (uint8_t *)realloc(out, (size_t)cap);
            prev = (uint8_t *)realloc(prev, (size_t)cap);
        }
        memcpy(out + noutb, entry, (size_t)elen);
        noutb += elen;
        if (nextCode <= 4095) {
            dict[nextCode] = (uint8_t *)ewcr_xmalloc((size_t)prevn + 1);
            memcpy(dict[nextCode], prev, (size_t)prevn);
            dict[nextCode][prevn] = entry[0];
            dlen[nextCode] = prevn + 1;
            nextCode++;
        }
        memcpy(prev, entry, (size_t)elen);
        prevn = elen;
        (void)tmp;
    }
    for (int k = 0; k < 4096; k++) free(dict[k]);
    free(dict);
    free(dlen);
    free(codes);
    free(prev);
    *nout = noutb;
    return out;
}

typedef struct {
    uint16_t *q;
    int64_t *sum;
    uint32_t *cnt;
    int n;
    int Q, Nframes;
} McMat;

static McMat spectral_variation(const int64_t *A, int nFrames, int Q, const double *gamma)
{
    int64_t *MAXV = (int64_t *)ewcr_xmalloc((size_t)Q * sizeof(int64_t));
    int64_t *MINV = (int64_t *)ewcr_xmalloc((size_t)Q * sizeof(int64_t));
    int64_t *SUMV = (int64_t *)ewcr_xmalloc((size_t)Q * sizeof(int64_t));
    uint32_t *CNTV = (uint32_t *)ewcr_xmalloc((size_t)Q * sizeof(uint32_t));
    memcpy(MAXV, A, (size_t)Q * sizeof(int64_t));
    memcpy(MINV, A, (size_t)Q * sizeof(int64_t));
    memcpy(SUMV, A, (size_t)Q * sizeof(int64_t));
    for (int i = 0; i < Q; i++) CNTV[i] = 1;

    int cap = Q * 4 + 16;
    uint16_t *qOut = (uint16_t *)ewcr_xcalloc((size_t)cap, sizeof(uint16_t));
    int64_t *sumOut = (int64_t *)ewcr_xcalloc((size_t)cap, sizeof(int64_t));
    uint32_t *cntOut = (uint32_t *)ewcr_xcalloc((size_t)cap, sizeof(uint32_t));
    int nOut = 0;

    for (int i = 1; i < nFrames; i++) {
        const int64_t *NEW = A + (size_t)i * Q;
        int *nov = (int *)ewcr_xcalloc((size_t)Q, sizeof(int));
        int nNew = 0;
        for (int q = 0; q < Q; q++) {
            int64_t nmax = MAXV[q] > NEW[q] ? MAXV[q] : NEW[q];
            int64_t nmin = MINV[q] < NEW[q] ? MINV[q] : NEW[q];
            int64_t ad = nmax - nmin;
            if (ad < 0) ad = -ad;
            if ((double)ad > gamma[q]) {
                nov[q] = 1;
                nNew++;
            }
        }
        if (nOut + nNew + Q + 8 > cap) {
            cap = (nOut + nNew + Q) * 2 + 64;
            qOut = (uint16_t *)realloc(qOut, (size_t)cap * sizeof(uint16_t));
            sumOut = (int64_t *)realloc(sumOut, (size_t)cap * sizeof(int64_t));
            cntOut = (uint32_t *)realloc(cntOut, (size_t)cap * sizeof(uint32_t));
        }
        for (int q = 0; q < Q; q++) {
            if (nov[q]) {
                qOut[nOut] = (uint16_t)(q + 1); /* 1-based q like MATLAB */
                sumOut[nOut] = SUMV[q];
                cntOut[nOut] = CNTV[q];
                nOut++;
                MAXV[q] = MINV[q] = SUMV[q] = NEW[q];
                CNTV[q] = 1;
            } else {
                int64_t nmax = MAXV[q] > NEW[q] ? MAXV[q] : NEW[q];
                int64_t nmin = MINV[q] < NEW[q] ? MINV[q] : NEW[q];
                MAXV[q] = nmax;
                MINV[q] = nmin;
                SUMV[q] += NEW[q];
                CNTV[q] += 1;
            }
        }
        free(nov);
    }
    if (nOut + Q > cap) {
        cap = nOut + Q + 8;
        qOut = (uint16_t *)realloc(qOut, (size_t)cap * sizeof(uint16_t));
        sumOut = (int64_t *)realloc(sumOut, (size_t)cap * sizeof(int64_t));
        cntOut = (uint32_t *)realloc(cntOut, (size_t)cap * sizeof(uint32_t));
    }
    for (int q = 0; q < Q; q++) {
        qOut[nOut] = (uint16_t)(q + 1);
        sumOut[nOut] = SUMV[q];
        cntOut[nOut] = CNTV[q];
        nOut++;
    }
    free(MAXV);
    free(MINV);
    free(SUMV);
    free(CNTV);
    McMat mc;
    mc.q = qOut;
    mc.sum = sumOut;
    mc.cnt = cntOut;
    mc.n = nOut;
    mc.Q = Q;
    mc.Nframes = nFrames;
    return mc;
}

static uint8_t *serialize_mc(const McMat *mc, int *nbytes)
{
    int n = mc->n;
    uint8_t *bytes = (uint8_t *)ewcr_xmalloc((size_t)n * 14);
    int p = 0;
    for (int k = 0; k < n; k++) {
        uint16_t q = mc->q[k];
        bytes[p++] = (uint8_t)(q & 0xFF);
        bytes[p++] = (uint8_t)((q >> 8) & 0xFF);
        uint64_t s;
        memcpy(&s, &mc->sum[k], 8);
        for (int b = 0; b < 8; b++) bytes[p++] = (uint8_t)((s >> (8 * b)) & 0xFF);
        uint32_t c = mc->cnt[k];
        for (int b = 0; b < 4; b++) bytes[p++] = (uint8_t)((c >> (8 * b)) & 0xFF);
    }
    *nbytes = n * 14;
    return bytes;
}

static McMat deserialize_mc(const uint8_t *bytes, int nbytes, int Q, int Nframes)
{
    int n = nbytes / 14;
    McMat mc;
    mc.q = (uint16_t *)ewcr_xmalloc((size_t)n * sizeof(uint16_t));
    mc.sum = (int64_t *)ewcr_xmalloc((size_t)n * sizeof(int64_t));
    mc.cnt = (uint32_t *)ewcr_xmalloc((size_t)n * sizeof(uint32_t));
    mc.n = n;
    mc.Q = Q;
    mc.Nframes = Nframes;
    int p = 0;
    for (int k = 0; k < n; k++) {
        mc.q[k] = (uint16_t)(bytes[p] | (bytes[p + 1] << 8));
        p += 2;
        uint64_t s = 0;
        for (int b = 0; b < 8; b++) s |= ((uint64_t)bytes[p++]) << (8 * b);
        memcpy(&mc.sum[k], &s, 8);
        uint32_t c = 0;
        for (int b = 0; b < 4; b++) c |= ((uint32_t)bytes[p++]) << (8 * b);
        mc.cnt[k] = c;
    }
    return mc;
}

static void reconstruct_MR(const McMat *mc, double *MR)
{
    int Q = mc->Q, Nf = mc->Nframes;
    memset(MR, 0, (size_t)Nf * (size_t)Q * sizeof(double));
    int *pos = (int *)ewcr_xmalloc((size_t)Q * sizeof(int));
    for (int i = 0; i < Q; i++) pos[i] = 0;
    for (int e = 0; e < mc->n; e++) {
        int q = (int)mc->q[e] - 1;
        int c = (int)mc->cnt[e];
        if (q < 0 || q >= Q || c < 1) continue;
        int a = pos[q];
        int b = a + c;
        if (b > Nf) b = Nf;
        double mean = (double)mc->sum[e] / (double)c;
        for (int i = a; i < b; i++) MR[i * Q + q] = mean;
        pos[q] = b;
    }
    free(pos);
}

static void free_mc(McMat *mc)
{
    free(mc->q);
    free(mc->sum);
    free(mc->cnt);
    mc->q = NULL;
}

int svdcs_codec_run(const double *x, int n, const SvdcsOpts *opt, double *xhat, EwcrRunResult *out)
{
    memset(out, 0, sizeof(*out));
    strcpy(out->algo, "SVDCS");
    SvdcsOpts D;
    memset(&D, 0, sizeof(D));
    D.fs = 12800;
    D.f0 = 50;
    D.Nppc = 256;
    D.signal_scale = 1;
    D.G = 0.03;
    D.beta = 0.10;
    D.input_full_scale_pu = 2;
    D.fft_full_scale_pu = 2;
    D.lagrange_order = 3;
    D.input_quantize = 1;
    D.input_bits = 16;
    D.fft_quantize = 1;
    D.fft_bits = 16;
    D.original_bits_per_sample = 16;
    D.use_lzw = 1;
    strcpy(D.sync_mode, "nominal");
    if (opt) D = *opt;

    double t0 = ewcr_now_s();
    double *xpu = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) xpu[i] = x[i] / D.signal_scale;
    if (D.input_quantize) {
        int32_t *qc = (int32_t *)ewcr_xmalloc((size_t)n * sizeof(int32_t));
        double step;
        quantize_signed(xpu, n, D.input_bits, D.input_full_scale_pu, qc, &step);
        for (int i = 0; i < n; i++) xpu[i] = (double)qc[i] * step;
        free(qc);
    }
    int M = D.Nppc;
    int maxF = n / (M > 0 ? (M / 2 > 0 ? 1 : 1) : 1) + 8;
    if (maxF < 4) maxF = 4;
    double *frames = (double *)ewcr_xmalloc((size_t)maxF * (size_t)M * sizeof(double));
    int nF = 0;
    if (prepare_cycles_nominal(xpu, n, &D, frames, &nF) != 0) {
        snprintf(out->note, sizeof(out->note), "sync failed");
        free(xpu);
        free(frames);
        return -1;
    }
    int Q = 2 * M;
    double *Afloat = (double *)ewcr_xcalloc((size_t)nF * Q, sizeof(double));
    cpx *row = (cpx *)ewcr_xmalloc((size_t)M * sizeof(cpx));
    for (int k = 0; k < nF; k++) {
        for (int i = 0; i < M; i++) row[i] = cpx_make(frames[k * M + i], 0);
        ewcr_fft(row, M, 0);
        for (int i = 0; i < M; i++) {
            cpx z = cpx_scale(row[i], 1.0 / M);
            Afloat[k * Q + 2 * i] = z.re;
            Afloat[k * Q + 2 * i + 1] = z.im;
        }
    }
    int32_t *Acode = (int32_t *)ewcr_xmalloc((size_t)nF * Q * sizeof(int32_t));
    double fftStep;
    if (D.fft_quantize) {
        quantize_signed(Afloat, nF * Q, D.fft_bits, D.fft_full_scale_pu, Acode, &fftStep);
    } else {
        fftStep = 1e-12;
        for (int i = 0; i < nF * Q; i++) Acode[i] = (int32_t)round(Afloat[i] / fftStep);
    }
    double *gamma = (double *)ewcr_xmalloc((size_t)Q * sizeof(double));
    for (int q = 0; q < Q; q++) {
        double g = D.G * pow((double)(q + 1), -D.beta);
        gamma[q] = g / fftStep;
    }
    int64_t *A64 = (int64_t *)ewcr_xmalloc((size_t)nF * Q * sizeof(int64_t));
    for (int i = 0; i < nF * Q; i++) A64[i] = (int64_t)Acode[i];
    McMat mc = spectral_variation(A64, nF, Q, gamma);
    int raw_n = 0;
    uint8_t *mcRaw = serialize_mc(&mc, &raw_n);
    int payload_n = 0, nCodes = 0;
    uint8_t *payload;
    if (D.use_lzw) {
        payload = lzw_encode12(mcRaw, raw_n, &payload_n, &nCodes);
        int chk_n = 0;
        uint8_t *chk = lzw_decode12(payload, nCodes, &chk_n);
        if (chk_n != raw_n || memcmp(chk, mcRaw, (size_t)raw_n) != 0) {
            snprintf(out->note, sizeof(out->note), "LZW roundtrip failed");
            out->ok = 0;
        }
        free(chk);
    } else {
        payload = mcRaw;
        payload_n = raw_n;
        mcRaw = NULL;
    }
    double bits_total = 8.0 * fmax(1.0, (double)payload_n);
    double t1 = ewcr_now_s();

    /* decode from payload */
    int raw2_n = 0;
    uint8_t *raw2;
    if (D.use_lzw) raw2 = lzw_decode12(payload, nCodes, &raw2_n);
    else {
        raw2 = payload;
        raw2_n = payload_n;
        payload = NULL;
    }
    McMat mc2 = deserialize_mc(raw2, raw2_n, Q, nF);
    double *MR = (double *)ewcr_xmalloc((size_t)nF * Q * sizeof(double));
    reconstruct_MR(&mc2, MR);
    int nsync = nF * M;
    double *xsync = (double *)ewcr_xmalloc((size_t)nsync * sizeof(double));
    for (int k = 0; k < nF; k++) {
        for (int i = 0; i < M; i++) row[i] = cpx_make(MR[k * Q + 2 * i] * fftStep, MR[k * Q + 2 * i + 1] * fftStep);
        for (int i = 0; i < M; i++) row[i] = cpx_scale(row[i], (double)M);
        ewcr_fft(row, M, 1);
        for (int i = 0; i < M; i++) xsync[k * M + i] = row[i].re * D.signal_scale;
    }
    /* compare against synchronized reference */
    int ncmp = nsync < n ? nsync : n;
    /* pad/trim xhat to original n: place sync samples, rest 0 */
    memset(xhat, 0, (size_t)n * sizeof(double));
    if (ncmp > n) ncmp = n;
    /* For metrics, MATLAB compares xref (sync frames) vs xhat (same length).
       We copy reconstructed sync into xhat[0:nsync] and will metric vs generated
       sync reference stored in frames. */
    memcpy(xhat, xsync, (size_t)ncmp * sizeof(double));
    double t2 = ewcr_now_s();

    double *xref = (double *)ewcr_xmalloc((size_t)nsync * sizeof(double));
    for (int i = 0; i < nsync; i++) xref[i] = frames[i] * D.signal_scale;
    /* If original n differs, report on sync domain (paper). For unified table use min length vs x. */
    if (nsync == n) {
        ewcr_unified_metrics(x, xhat, n, bits_total, D.fs, D.f0, D.original_bits_per_sample, &out->metrics);
    } else {
        /* metric vs synchronized reference, but CR vs original n*bits? MATLAB uses xref length. */
        ewcr_unified_metrics(xref, xsync, nsync, bits_total, D.fs, D.f0, D.original_bits_per_sample, &out->metrics);
        memcpy(xhat, xsync, (size_t)((nsync < n ? nsync : n) * sizeof(double)));
    }
    out->enc_s = t1 - t0;
    out->dec_s = t2 - t1;
    out->N = out->metrics.N;
    if (!out->ok && out->note[0]) {
        /* LZW fail */
    } else {
        out->ok = 1;
        snprintf(out->note, sizeof(out->note), "nF=%d Nppc=%d payload=%d LZW=%d", nF, M, payload_n, nCodes);
    }

    free(xpu);
    free(frames);
    free(Afloat);
    free(row);
    free(Acode);
    free(gamma);
    free(A64);
    free_mc(&mc);
    free(mcRaw);
    free(payload);
    if (D.use_lzw) free(raw2);
    free_mc(&mc2);
    free(MR);
    free(xsync);
    free(xref);
    return out->ok ? 0 : -1;
}
