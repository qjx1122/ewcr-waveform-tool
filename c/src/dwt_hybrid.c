#include "ewcr.h"

#include <stdlib.h>
#include <string.h>

static double soft_thresh(double v, double T)
{
    double a = fabs(v) - T;
    if (a < 0) a = 0;
    return (v >= 0 ? 1.0 : -1.0) * a;
}

static void energy_retain(double *bj, int n, double frac)
{
    if (n <= 0) return;
    int *ord = (int *)ewcr_xmalloc((size_t)n * sizeof(int));
    double *ab = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
    double tot = 0.0;
    for (int i = 0; i < n; i++) {
        ord[i] = i;
        ab[i] = fabs(bj[i]);
        tot += bj[i] * bj[i];
    }
    if (tot == 0.0) {
        free(ord);
        free(ab);
        return;
    }
    for (int i = 0; i < n; i++) {
        int bi = i;
        for (int j = i + 1; j < n; j++) if (ab[ord[j]] > ab[ord[bi]]) bi = j;
        int t = ord[i];
        ord[i] = ord[bi];
        ord[bi] = t;
    }
    double cum = 0.0, need = frac * tot;
    int kkeep = n;
    for (int i = 0; i < n; i++) {
        double v = bj[ord[i]];
        cum += v * v;
        if (cum >= need) {
            kkeep = i + 1;
            break;
        }
    }
    int *keep = (int *)ewcr_xcalloc((size_t)n, sizeof(int));
    for (int i = 0; i < kkeep; i++) keep[ord[i]] = 1;
    for (int i = 0; i < n; i++) if (!keep[i]) bj[i] = 0.0;
    free(ord);
    free(ab);
    free(keep);
}

static double median_abs(const double *bj, int n)
{
    double *t = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) t[i] = fabs(bj[i]);
    for (int i = 1; i < n; i++) {
        double v = t[i];
        int j = i - 1;
        while (j >= 0 && t[j] > v) {
            t[j + 1] = t[j];
            j--;
        }
        t[j + 1] = v;
    }
    double m = (n % 2) ? t[n / 2] : 0.5 * (t[n / 2 - 1] + t[n / 2]);
    free(t);
    return m;
}

static double quant_step(const double *blk, int n, int qbits)
{
    double mx = 0.0;
    for (int i = 0; i < n; i++) {
        double a = fabs(blk[i]);
        if (a > mx) mx = a;
    }
    if (mx == 0.0) return 1.0;
    return mx / (double)(1 << (qbits - 1));
}

static void quantize_block(const double *blk, int n, double qs, int qbits, double *vq)
{
    double lo = -(double)(1 << (qbits - 1));
    double hi = (double)(1 << (qbits - 1)) - 1.0;
    for (int i = 0; i < n; i++) {
        double v = round(blk[i] / qs);
        if (v < lo) v = lo;
        if (v > hi) v = hi;
        vq[i] = v;
    }
}

int dwt_codec_run(const double *x, int n, const DwtOpts *opt, double *xhat, EwcrRunResult *out)
{
    memset(out, 0, sizeof(*out));
    strcpy(out->algo, "DWT-Hybrid");
    DwtOpts D;
    memset(&D, 0, sizeof(D));
    strcpy(D.thr_mode, "energy");
    D.thr_frac = 0.99;
    D.qbits = 8;
    D.delta = 1;
    D.RQ = 16;
    if (opt) D = *opt;

    double t0 = ewcr_now_s();
    int L;
    if (D.level <= 0) L = (int)floor(log2(n / 16.0));
    else L = D.level;
    if (L > 10) L = 10;
    if (L < 1) L = 1;
    while (L > 1 && (n >> L) * (1 << L) != n) L--;
    if ((n >> L) << L != n) {
        /* trim L until n divisible by 2^L */
        while (L > 1 && (n % (1 << L)) != 0) L--;
    }

    double *C = (double *)ewcr_xcalloc((size_t)n, sizeof(double));
    int book[32];
    int ncoef = 0;
    if (ewcr_wavedec_db4(x, n, L, C, book, &ncoef) != 0) {
        snprintf(out->note, sizeof(out->note), "wavedec failed N=%d L=%d", n, L);
        free(C);
        return -1;
    }
    int nA = book[0];
    double *Cq = (double *)ewcr_xcalloc((size_t)ncoef, sizeof(double));
    double *scales = (double *)ewcr_xcalloc((size_t)(L + 1), sizeof(double));

    scales[0] = quant_step(C, nA, D.qbits);
    quantize_block(C, nA, scales[0], D.qbits, Cq);

    int pos = nA;
    int univ = (strcmp(D.thr_mode, "universal") == 0);
    for (int j = 0; j < L; j++) {
        int len = book[j + 1];
        double *bj = (double *)ewcr_xmalloc((size_t)len * sizeof(double));
        memcpy(bj, C + pos, (size_t)len * sizeof(double));
        if (univ) {
            double sigma = median_abs(bj, len) / 0.6745;
            double Tj = sigma * sqrt(2.0 * log((double)len));
            for (int i = 0; i < len; i++) bj[i] = soft_thresh(bj[i], Tj);
        } else {
            energy_retain(bj, len, D.thr_frac);
        }
        scales[j + 1] = quant_step(bj, len, D.qbits);
        quantize_block(bj, len, scales[j + 1], D.qbits, Cq + pos);
        pos += len;
        free(bj);
    }

    double *dA = (double *)ewcr_xmalloc((size_t)nA * sizeof(double));
    if (D.delta) {
        dA[0] = Cq[0];
        for (int i = 1; i < nA; i++) dA[i] = Cq[i] - Cq[i - 1];
    } else {
        memcpy(dA, Cq, (size_t)nA * sizeof(double));
    }

    int nD = ncoef - nA;
    int32_t *allvals = (int32_t *)ewcr_xmalloc((size_t)nD * sizeof(int32_t));
    for (int i = 0; i < nD; i++) allvals[i] = (int32_t)Cq[nA + i];

    int nz_count = 0;
    for (int i = 0; i < nD; i++) if (allvals[i] != 0) nz_count++;
    int ngaps = nz_count + 1;
    int32_t *gaps = (int32_t *)ewcr_xmalloc((size_t)ngaps * sizeof(int32_t));
    if (nz_count == 0) {
        gaps[0] = nD;
    } else {
        int gi = 0;
        int last = -1;
        for (int i = 0; i < nD; i++) {
            if (allvals[i] != 0) {
                gaps[gi++] = i - last - 1;
                last = i;
            }
        }
        gaps[gi] = nD - 1 - last;
    }

    int32_t *dAi = (int32_t *)ewcr_xmalloc((size_t)nA * sizeof(int32_t));
    for (int i = 0; i < nA; i++) dAi[i] = (int32_t)dA[i];

    uint8_t *b_vals = NULL, *cb_v = NULL, *b_gaps = NULL, *cb_g = NULL, *b_a = NULL, *cb_a = NULL;
    int32_t *sy_v = NULL, *sy_g = NULL, *sy_a = NULL;
    int *cl_v = NULL, *cl_g = NULL, *cl_a = NULL;
    int nb_v = 0, nb_g = 0, nb_a = 0, ns_v = 0, ns_g = 0, ns_a = 0;
    ewcr_huff_encode(allvals, nD, &b_vals, &nb_v, &sy_v, &ns_v, &cb_v, &cl_v);
    ewcr_huff_encode(gaps, ngaps, &b_gaps, &nb_g, &sy_g, &ns_g, &cb_g, &cl_g);
    ewcr_huff_encode(dAi, nA, &b_a, &nb_a, &sy_a, &ns_a, &cb_a, &cl_a);

    /* MATLAB: nbits = numel(enco) + numel(syms)*48 */
    double bits_vals = nb_v + ns_v * 48.0;
    double bits_gaps = nb_g + ns_g * 48.0;
    double bits_a = nb_a + ns_a * 48.0;
    double bits_header = 128 + (L + 1) * 32 + 64;
    double bits_total = bits_vals + bits_gaps + bits_a + bits_header;
    double t1 = ewcr_now_s();

    /* Huffman decode path */
    int32_t *dA_dec = (int32_t *)ewcr_xcalloc((size_t)nA, sizeof(int32_t));
    int32_t *vals_dec = (int32_t *)ewcr_xcalloc((size_t)nD, sizeof(int32_t));
    int32_t *gaps_dec = (int32_t *)ewcr_xcalloc((size_t)ngaps, sizeof(int32_t));
    ewcr_huff_decode(b_a, nb_a, sy_a, ns_a, cb_a, cl_a, dA_dec, nA);
    ewcr_huff_decode(b_vals, nb_v, sy_v, ns_v, cb_v, cl_v, vals_dec, nD);
    ewcr_huff_decode(b_gaps, nb_g, sy_g, ns_g, cb_g, cl_g, gaps_dec, ngaps);

    double *Aq = (double *)ewcr_xmalloc((size_t)nA * sizeof(double));
    if (D.delta) {
        Aq[0] = dA_dec[0];
        for (int i = 1; i < nA; i++) Aq[i] = Aq[i - 1] + dA_dec[i];
    } else {
        for (int i = 0; i < nA; i++) Aq[i] = dA_dec[i];
    }

    double *Cq2 = (double *)ewcr_xcalloc((size_t)ncoef, sizeof(double));
    for (int i = 0; i < nA; i++) Cq2[i] = Aq[i];
    /* nonzero values from decoded allvals, placed by gaps (MATLAB decoder) */
    int nvals = 0;
    for (int i = 0; i < nD; i++) if (vals_dec[i] != 0) nvals++;
    int32_t *nzv = (int32_t *)ewcr_xmalloc((size_t)(nvals + 1) * sizeof(int32_t));
    nvals = 0;
    for (int i = 0; i < nD; i++) if (vals_dec[i] != 0) nzv[nvals++] = vals_dec[i];
    pos = nA; /* 0-based start of details; MATLAB pos = nA+1 then pos += gaps; place; pos++ */
    /* MATLAB: pos = nA+1 (1-based) => 0-based nA
       pos = pos + gaps(gi); place at pos; pos = pos+1
       In 0-based: pos starts nA, then pos += gaps, write Cq[pos], pos++ */
    int gi = 0;
    int p = nA;
    for (int i = 0; i < nvals; i++) {
        if (gi >= ngaps) break;
        p += gaps_dec[gi];
        gi++;
        if (p >= ncoef) break;
        Cq2[p] = (double)nzv[i];
        p += 1;
    }

    double *Crec = (double *)ewcr_xcalloc((size_t)ncoef, sizeof(double));
    for (int i = 0; i < nA; i++) Crec[i] = Cq2[i] * scales[0];
    pos = nA;
    for (int j = 0; j < L; j++) {
        int len = book[j + 1];
        for (int i = 0; i < len; i++) Crec[pos + i] = Cq2[pos + i] * scales[j + 1];
        pos += len;
    }
    ewcr_waverec_db4(Crec, book, L, xhat);
    double t2 = ewcr_now_s();

    double fs = (D.fs > 0.0) ? D.fs : 12800.0;
    double f0 = (D.f0 > 0.0) ? D.f0 : 50.0;
    ewcr_unified_metrics(x, xhat, n, bits_total, fs, f0, D.RQ, &out->metrics);
    out->enc_s = t1 - t0;
    out->dec_s = t2 - t1;
    out->N = n;
    out->ok = 1;
    snprintf(out->note, sizeof(out->note), "L=%d nA=%d nz=%d", L, nA, nz_count);

    {
        EwcrBuf b;
        ewcr_buf_init(&b);
        ewcr_pack_header(&b, 2, n, fs, f0, bits_total);
        ewcr_buf_u32(&b, (uint32_t)L);
        ewcr_buf_u32(&b, (uint32_t)D.qbits);
        ewcr_buf_u32(&b, (uint32_t)D.delta);
        ewcr_buf_u32(&b, (uint32_t)nA);
        ewcr_buf_u32(&b, (uint32_t)ncoef);
        for (int i = 0; i <= L + 1; i++) ewcr_buf_i32(&b, book[i]);
        for (int i = 0; i <= L; i++) ewcr_buf_f64(&b, scales[i]);
        for (int i = 0; i < ncoef; i++) ewcr_buf_i32(&b, (int32_t)Cq[i]);
        out->compressed = b.d;
        out->compressed_nbytes = (int)b.n;
    }

    free(C);
    free(Cq);
    free(scales);
    free(dA);
    free(allvals);
    free(gaps);
    free(dAi);
    free(b_vals);
    free(cb_v);
    free(sy_v);
    free(cl_v);
    free(b_gaps);
    free(cb_g);
    free(sy_g);
    free(cl_g);
    free(b_a);
    free(cb_a);
    free(sy_a);
    free(cl_a);
    free(dA_dec);
    free(vals_dec);
    free(gaps_dec);
    free(Aq);
    free(Cq2);
    free(nzv);
    free(Crec);
    return 0;
}
