#include "ewcr.h"

#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    int *idx;
    int nidx;
    double *Re, *Im;
    double qe;
} AsbcEfft;

static cpx quant_uni(cpx v, int R, double qstep)
{
    double cr = round((v.re + 1.0) / qstep);
    double ci = round((v.im + 1.0) / qstep);
    double mx = (double)((1 << R) - 1);
    if (cr < 0) cr = 0;
    if (cr > mx) cr = mx;
    if (ci < 0) ci = 0;
    if (ci > mx) ci = mx;
    return cpx_make(cr, ci);
}

static cpx dequant_uni(cpx c, double qstep)
{
    return cpx_make(c.re * qstep - 1.0, c.im * qstep - 1.0);
}

int asbc_codec_run(const double *x, int n, double fs, const AsbcOpts *opt,
                   double *xhat, EwcrRunResult *out)
{
    memset(out, 0, sizeof(*out));
    strcpy(out->algo, "ASBC");
    AsbcOpts D = {20, 1e-4, 0.99, 8, 5, 48, 16};
    if (opt) D = *opt;
    double t0 = ewcr_now_s();

    int fft_len = 1 << ewcr_nextpow2(n);
    if (fft_len < n) fft_len = n;
    cpx *Xf = (cpx *)ewcr_xmalloc((size_t)fft_len * sizeof(cpx));
    ewcr_fft_real(x, n, fft_len, Xf);

    double *mag = (double *)ewcr_xcalloc((size_t)fft_len, sizeof(double));
    int im = 0;
    double best = -1.0;
    for (int i = 0; i < fft_len; i++) {
        mag[i] = sqrt(cpx_abs2(Xf[i]));
        double freq = (double)i / fft_len * fs;
        if (freq >= 45.0 && freq <= 55.0 && mag[i] > best) {
            best = mag[i];
            im = i;
        }
    }
    double F0 = (double)im / fft_len * fs;
    if (F0 <= 0) F0 = 50.0;
    double w0 = 2.0 * M_PI * F0 / fs;

    int blk = (int)round(D.block_cycles * fs / F0);
    int nb = (blk > 0) ? (n / blk) : 1;
    if (nb < 2) {
        blk = n;
        nb = 1;
    }

    double W = D.W;
    int Kmax = (int)floor((fs / 2.0 - W / 2.0) / F0);
    if (Kmax < 1) Kmax = 1;
    double bin_spacing = fs / (double)fft_len;
    double *Ek = (double *)ewcr_xcalloc((size_t)(Kmax + 1), sizeof(double));
    for (int k = 1; k <= Kmax; k++) {
        int kc = (int)round(k * F0 / bin_spacing); /* 0-based = MATLAB round + 1 - 1 */
        if (kc >= 0 && kc < fft_len) Ek[k] = mag[kc] * mag[kc];
    }
    double thrE = D.act_thresh * Ek[1];
    int K = 1;
    for (int k = 2; k <= Kmax; k++) if (Ek[k] > thrE) K = k;

    double scale = 0.0;
    for (int i = 0; i < n; i++) {
        double a = fabs(x[i]);
        if (a > scale) scale = a;
    }
    if (scale == 0.0) scale = 1.0;
    double *xnorm = (double *)ewcr_xmalloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) xnorm[i] = x[i] / scale;
    double qstep = 2.0 / (double)(1 << D.Rk);
    int Sk = (int)ceil(fs / W);
    if (Sk < 1) Sk = 1;

    int *keep = (int *)ewcr_xcalloc((size_t)blk, sizeof(int));
    for (int i = 0; i < blk; i++) {
        if (i <= (W / 2.0) / fs * blk || i >= blk - (W / 2.0) / fs * blk) keep[i] = 1;
    }

    int *mask_h = (int *)ewcr_xcalloc((size_t)nb * K, sizeof(int));
    int *mask_e = (int *)ewcr_xcalloc((size_t)nb, sizeof(int));
    int *sub_len = (int *)ewcr_xcalloc((size_t)nb * K, sizeof(int));
    cpx **sub_data = (cpx **)ewcr_xcalloc((size_t)nb * K, sizeof(cpx *));
    AsbcEfft *e_fft = (AsbcEfft *)ewcr_xcalloc((size_t)nb, sizeof(AsbcEfft));

    cpx *yb = (cpx *)ewcr_xmalloc((size_t)blk * sizeof(cpx));
    cpx *ylp = (cpx *)ewcr_xmalloc((size_t)blk * sizeof(cpx));
    cpx *yup = (cpx *)ewcr_xmalloc((size_t)blk * sizeof(cpx));
    double *xrec = (double *)ewcr_xmalloc((size_t)blk * sizeof(double));
    double *eb = (double *)ewcr_xmalloc((size_t)blk * sizeof(double));

    for (int b = 0; b < nb; b++) {
        const double *xb = xnorm + b * blk;
        /* fundamental always active */
        for (int i = 0; i < blk; i++) {
            double ang = -w0 * i;
            cpx e = cpx_make(cos(ang), sin(ang));
            yb[i] = cpx_mul(cpx_make(xb[i], 0), e);
        }
        memcpy(ylp, yb, (size_t)blk * sizeof(cpx));
        ewcr_fft(ylp, blk, 0);
        for (int i = 0; i < blk; i++) if (!keep[i]) ylp[i] = cpx_make(0, 0);
        ewcr_fft(ylp, blk, 1);
        int nds = 0;
        for (int i = 0; i < blk; i += Sk) nds++;
        cpx *yds = (cpx *)ewcr_xmalloc((size_t)nds * sizeof(cpx));
        int t = 0;
        for (int i = 0; i < blk; i += Sk) yds[t++] = quant_uni(ylp[i], D.Rk, qstep);
        sub_data[b * K + 0] = yds;
        sub_len[b * K + 0] = nds;
        mask_h[b * K + 0] = 1;

        double exb = 0.0;
        for (int i = 0; i < blk; i++) exb += xb[i] * xb[i];

        for (int k = 2; k <= K; k++) {
            for (int i = 0; i < blk; i++) {
                double ang = -k * w0 * i;
                cpx e = cpx_make(cos(ang), sin(ang));
                yb[i] = cpx_mul(cpx_make(xb[i], 0), e);
            }
            memcpy(ylp, yb, (size_t)blk * sizeof(cpx));
            ewcr_fft(ylp, blk, 0);
            for (int i = 0; i < blk; i++) if (!keep[i]) ylp[i] = cpx_make(0, 0);
            ewcr_fft(ylp, blk, 1);
            double ey = 0.0;
            for (int i = 0; i < blk; i++) ey += cpx_abs2(ylp[i]);
            if (ey > D.act_thresh * exb) {
                nds = 0;
                for (int i = 0; i < blk; i += Sk) nds++;
                yds = (cpx *)ewcr_xmalloc((size_t)nds * sizeof(cpx));
                t = 0;
                for (int i = 0; i < blk; i += Sk) yds[t++] = quant_uni(ylp[i], D.Rk, qstep);
                sub_data[b * K + (k - 1)] = yds;
                sub_len[b * K + (k - 1)] = nds;
                mask_h[b * K + (k - 1)] = 1;
            }
        }

        memset(xrec, 0, (size_t)blk * sizeof(double));
        for (int k = 1; k <= K; k++) {
            if (!mask_h[b * K + (k - 1)]) continue;
            cpx *qd = sub_data[b * K + (k - 1)];
            int ns = sub_len[b * K + (k - 1)];
            cpx *deq = (cpx *)ewcr_xmalloc((size_t)ns * sizeof(cpx));
            for (int i = 0; i < ns; i++) deq[i] = dequant_uni(qd[i], qstep);
            ewcr_interpft(deq, ns, yup, blk);
            free(deq);
            for (int i = 0; i < blk; i++) {
                double ang = k * w0 * i;
                cpx e = cpx_make(cos(ang), sin(ang));
                cpx p = cpx_mul(yup[i], e);
                xrec[i] += sqrt(2.0) * p.re;
            }
        }
        for (int i = 0; i < blk; i++) eb[i] = xb[i] - xrec[i];
        double eeb = 0.0;
        for (int i = 0; i < blk; i++) eeb += eb[i] * eb[i];
        if (eeb > D.act_thresh * exb) {
            mask_e[b] = 1;
            cpx *Eb = (cpx *)ewcr_xcalloc((size_t)blk, sizeof(cpx));
            for (int i = 0; i < blk; i++) Eb[i] = cpx_make(eb[i], 0);
            ewcr_fft(Eb, blk, 0);
            int half = blk / 2;
            int namp = half; /* bins 1..half  (0-based) */
            double *amp = (double *)ewcr_xmalloc((size_t)namp * sizeof(double));
            int *ord = (int *)ewcr_xmalloc((size_t)namp * sizeof(int));
            double sumamp2 = 0.0, maxamp = 0.0;
            for (int i = 0; i < namp; i++) {
                amp[i] = sqrt(cpx_abs2(Eb[i + 1]));
                ord[i] = i;
                sumamp2 += amp[i] * amp[i];
                if (amp[i] > maxamp) maxamp = amp[i];
            }
            /* argsort descend */
            for (int i = 0; i < namp; i++) {
                int bi = i;
                for (int j = i + 1; j < namp; j++) if (amp[ord[j]] > amp[ord[bi]]) bi = j;
                int tmp = ord[i];
                ord[i] = ord[bi];
                ord[bi] = tmp;
            }
            double cum = 0.0, need = D.energy_frac * sumamp2;
            int kk = namp;
            for (int i = 0; i < namp; i++) {
                cum += amp[ord[i]] * amp[ord[i]];
                if (cum >= need) {
                    kk = i + 1;
                    break;
                }
            }
            if (kk > D.ke_max) kk = D.ke_max;
            int *sel = (int *)ewcr_xmalloc((size_t)kk * sizeof(int));
            for (int i = 0; i < kk; i++) sel[i] = ord[i] + 1; /* MATLAB sel: 1-based amp index = 0-based bin */
            /* sort sel */
            for (int i = 1; i < kk; i++) {
                int v = sel[i], j = i - 1;
                while (j >= 0 && sel[j] > v) {
                    sel[j + 1] = sel[j];
                    j--;
                }
                sel[j + 1] = v;
            }
            double qe = maxamp / (double)(1 << 15);
            if (qe == 0.0) qe = 1.0;
            e_fft[b].nidx = kk;
            e_fft[b].idx = sel;
            e_fft[b].qe = qe;
            e_fft[b].Re = (double *)ewcr_xmalloc((size_t)kk * sizeof(double));
            e_fft[b].Im = (double *)ewcr_xmalloc((size_t)kk * sizeof(double));
            for (int i = 0; i < kk; i++) {
                int bin = sel[i]; /* 0-based bin */
                e_fft[b].Re[i] = round(Eb[bin].re / qe);
                e_fft[b].Im[i] = round(Eb[bin].im / qe);
            }
            free(Eb);
            free(amp);
            free(ord);
        }
    }

    /* bit accounting */
    double bits_sub = 0;
    for (int b = 0; b < nb; b++)
        for (int k = 0; k < K; k++)
            if (mask_h[b * K + k]) bits_sub += 2.0 * sub_len[b * K + k] * D.Rk;
    double bits_mask = (double)nb * K + nb;
    int idx_bits = (int)ceil(log2(blk / 2.0));
    if (idx_bits < 1) idx_bits = 1;
    double bits_ie = 0;
    for (int b = 0; b < nb; b++) if (mask_e[b]) {
        int kk = e_fft[b].nidx;
        bits_ie += 32 + kk * (idx_bits + 2 * 32);
    }
    double bits_total = 128 + 64 + bits_mask + bits_sub + bits_ie;
    double t1 = ewcr_now_s();

    /* decode */
    memset(xhat, 0, (size_t)n * sizeof(double));
    for (int b = 0; b < nb; b++) {
        memset(xrec, 0, (size_t)blk * sizeof(double));
        for (int k = 1; k <= K; k++) {
            if (!mask_h[b * K + (k - 1)]) continue;
            cpx *qd = sub_data[b * K + (k - 1)];
            int ns = sub_len[b * K + (k - 1)];
            cpx *deq = (cpx *)ewcr_xmalloc((size_t)ns * sizeof(cpx));
            for (int i = 0; i < ns; i++) deq[i] = dequant_uni(qd[i], qstep);
            ewcr_interpft(deq, ns, yup, blk);
            free(deq);
            for (int i = 0; i < blk; i++) {
                double ang = k * w0 * i;
                cpx e = cpx_make(cos(ang), sin(ang));
                cpx p = cpx_mul(yup[i], e);
                xrec[i] += sqrt(2.0) * p.re;
            }
        }
        if (mask_e[b]) {
            cpx *Eb = (cpx *)ewcr_xcalloc((size_t)blk, sizeof(cpx));
            for (int i = 0; i < e_fft[b].nidx; i++) {
                int sel = e_fft[b].idx[i]; /* 0-based bin */
                cpx v = cpx_make(e_fft[b].Re[i] * e_fft[b].qe, e_fft[b].Im[i] * e_fft[b].qe);
                Eb[sel] = v;
                if (sel > 0 && sel < blk) Eb[blk - sel] = cpx_conj(v);
            }
            ewcr_fft(Eb, blk, 1);
            for (int i = 0; i < blk; i++) xrec[i] += Eb[i].re;
            free(Eb);
        }
        for (int i = 0; i < blk; i++) xhat[b * blk + i] = xrec[i] * scale;
    }
    double t2 = ewcr_now_s();

    ewcr_unified_metrics(x, xhat, n, bits_total, fs, 50.0, D.RQ, &out->metrics);
    out->enc_s = t1 - t0;
    out->dec_s = t2 - t1;
    out->N = n;
    out->ok = 1;
    snprintf(out->note, sizeof(out->note), "F0=%.4f K=%d nb=%d blk=%d", F0, K, nb, blk);

    for (int i = 0; i < nb * K; i++) free(sub_data[i]);
    for (int b = 0; b < nb; b++) {
        free(e_fft[b].idx);
        free(e_fft[b].Re);
        free(e_fft[b].Im);
    }
    free(Xf);
    free(mag);
    free(Ek);
    free(xnorm);
    free(keep);
    free(mask_h);
    free(mask_e);
    free(sub_len);
    free(sub_data);
    free(e_fft);
    free(yb);
    free(ylp);
    free(yup);
    free(xrec);
    free(eb);
    return 0;
}
