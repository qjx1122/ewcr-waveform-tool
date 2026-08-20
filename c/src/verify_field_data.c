#include "ewcr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int file_exists(const char *p)
{
    FILE *f = fopen(p, "r");
    if (!f) return 0;
    fclose(f);
    return 1;
}

static const char *find_data_dir(void)
{
    static const char *cands[] = {"../data", "data", "../../data", "/home/user/ewcr-waveform-tool/data", NULL};
    static char buf[512];
    for (int i = 0; cands[i]; i++) {
        snprintf(buf, sizeof(buf), "%s/wave1.csv", cands[i]);
        if (file_exists(buf)) return cands[i];
    }
    return NULL;
}

static int finite_vec(const double *x, int n)
{
    for (int i = 0; i < n; i++) if (!isfinite(x[i])) return 0;
    return 1;
}

static double peak_abs(const double *x, int n)
{
    double m = 0.0;
    for (int i = 0; i < n; i++) {
        double a = fabs(x[i]);
        if (a > m) m = a;
    }
    return m;
}

static void apply_field_cfg(EwcrBenchCfg *cfg, double fs, double f0, double peak)
{
    cfg->fs = fs;
    cfg->f0 = f0;
    cfg->mmc.fs = fs;
    cfg->mmc.fn = f0;
    cfg->mmc.N = (int)round(fs / f0);
    if (cfg->mmc.N < 8) cfg->mmc.N = 8;
    cfg->mmc.n_tot = cfg->mmc.N;
    cfg->mmc.verbose = 0;
    cfg->svdcs.fs = fs;
    cfg->svdcs.f0 = f0;
    cfg->svdcs.Nppc = (int)round(fs / f0);
    cfg->svdcs.signal_scale = (peak > 1e-9) ? peak : 1.0;
    cfg->svdcs.verbose = 0;
}

int main(int argc, char **argv)
{
    int cycles = 10;
    if (argc > 1 && strcmp(argv[1], "--quick") == 0) cycles = 4;

    const char *datadir = find_data_dir();
    if (!datadir) {
        fprintf(stderr, "Cannot find data/wave1.csv (tried ../data, data, ...)\n");
        return 2;
    }
    printf("EWCR field-data five-codec verification\n");
    printf("data dir: %s\n", datadir);
    printf("segment: %d cycles @ estimated fs (skip 1 cycle)\n", cycles);
    printf("压缩比 CR = (N*16 bit)/bits_compressed;  相似度%% = 100*Pearson(original, reconstructed)\n\n");

    const char *files[] = {"wave1.csv", "wave2.csv", "wave3.csv", "wave4.csv"};
    const char *chans[] = {"UA", "IA"};
    const char *algos[] = {"ASBC", "DWT-Hybrid", "CS-OMP", "MMC", "SVDCS"};

    mkdir("results", 0755);
    mkdir("results/out", 0755);
    FILE *man = fopen("results/out/manifest.csv", "w");
    if (man)
        fprintf(man, "stem,file,channel,algo,original_csv,reconstructed_csv,compressed_ewcr,compressed_bytes\n");
    FILE *csv = fopen("results/verify_field_data.csv", "w");
    if (!csv) csv = fopen("c/results/verify_field_data.csv", "w");
    if (csv)
        fprintf(csv, "file,channel,algo,fs,N,ok,compress_ratio,similarity_pct,corr,SNR_dB,"
                     "NMSE_dB,RMSE,PRD,MAXE,PSNR_dB,enc_s,dec_s,note\n");

    printf("%-10s %-4s %-12s %7s %6s %10s %10s %9s %8s %8s  %s\n", "file", "ch", "algo", "fs", "N",
           "CR", "sim%", "SNR_dB", "enc_ms", "dec_ms", "status");
    printf("----------------------------------------------------------------------------------------------------\n");

    int nrun = 0, nfail = 0;
    for (int fi = 0; fi < 4; fi++) {
        char path[768];
        snprintf(path, sizeof(path), "%s/%s", datadir, files[fi]);
        EwcrWaveCsv w;
        if (ewcr_load_wave_csv(path, &w) != 0) {
            fprintf(stderr, "FAILED to load %s\n", path);
            nfail++;
            continue;
        }
        printf("# %s  n=%d  fs=%.2f Hz  dur=%.3f s  UA_peak=%.2f  IA_peak=%.3f\n", files[fi], w.n,
               w.fs, w.n / w.fs, peak_abs(w.UA, w.n), peak_abs(w.IA, w.n));

        double f0 = 50.0;
        int spc = (int)round(w.fs / f0);
        int skip = spc; /* skip first cycle */
        int nseg = cycles * spc;
        if (skip + nseg > w.n) {
            skip = 0;
            nseg = w.n;
            /* trim to whole cycles */
            nseg = (nseg / spc) * spc;
        }
        if (nseg < spc * 2) {
            fprintf(stderr, "  record too short for %d cycles\n", cycles);
            ewcr_free_wave_csv(&w);
            nfail++;
            continue;
        }

        double *xhat = (double *)ewcr_xmalloc((size_t)nseg * sizeof(double));
        for (int ci = 0; ci < 2; ci++) {
            const double *full = ewcr_wave_channel(&w, chans[ci]);
            const double *x = full + skip;
            double peak = peak_abs(x, nseg);
            EwcrBenchCfg cfg;
            ewcr_default_cfg(&cfg);
            apply_field_cfg(&cfg, w.fs, f0, peak);

            for (int a = 0; a < 5; a++) {
                memset(xhat, 0, (size_t)nseg * sizeof(double));
                EwcrRunResult r;
                int rc = ewcr_run_codec(algos[a], x, nseg, &cfg, xhat, &r);
                int nuse = r.N > 0 ? r.N : nseg;
                if (nuse > nseg) nuse = nseg;
                /* recompute metrics with field fs/f0 */
                if (rc == 0 && nuse > 0)
                    ewcr_unified_metrics(x, xhat, nuse, r.metrics.bits_compressed, w.fs, f0,
                                         cfg.raw_bits_per_sample, &r.metrics);
                int ok = (rc == 0 && r.ok && finite_vec(xhat, nuse) && isfinite(r.metrics.CR) &&
                          r.metrics.CR > 0 && isfinite(r.metrics.SNR_dB));
                r.ok = ok;
                nrun++;
                if (!ok) nfail++;
                printf("%-10s %-4s %-12s %7.0f %6d %10.3f %10.4f %9.2f %8.1f %8.1f  %s (%s)\n",
                       files[fi], chans[ci], algos[a], w.fs, r.N, r.metrics.CR,
                       r.metrics.similarity_pct, r.metrics.SNR_dB, 1000.0 * r.enc_s,
                       1000.0 * r.dec_s, ok ? "PASS" : "FAIL", r.note);
                fflush(stdout);
                if (csv)
                    fprintf(csv,
                            "%s,%s,%s,%.8g,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,\"%s\"\n",
                            files[fi], chans[ci], algos[a], w.fs, r.N, ok, r.metrics.CR,
                            r.metrics.similarity_pct, r.metrics.corr, r.metrics.SNR_dB,
                            r.metrics.NMSE_dB, r.metrics.RMSE, r.metrics.PRD, r.metrics.MAXE,
                            r.metrics.PSNR_dB, r.enc_s, r.dec_s, r.note);
                if (ok) {
                    char fileid[32];
                    snprintf(fileid, sizeof(fileid), "%s", files[fi]);
                    char *dot = strchr(fileid, '.');
                    if (dot) *dot = 0;
                    char algo_id[32];
                    snprintf(algo_id, sizeof(algo_id), "%s", algos[a]);
                    for (char *p = algo_id; *p; p++) if (*p == '-') *p = '_';
                    char stem[128];
                    snprintf(stem, sizeof(stem), "%s_%s_%s", fileid, chans[ci], algo_id);
                    ewcr_save_case_files("results/out", stem, x, xhat, nuse, w.fs, &r);
                    if (man)
                        fprintf(man, "%s,%s,%s,%s,results/out/%s_original.csv,results/out/%s_reconstructed.csv,results/out/%s_compressed.ewcr,%d\n",
                                stem, files[fi], chans[ci], algos[a], stem, stem, stem,
                                r.compressed_nbytes);
                }
                ewcr_result_release(&r);
            }
        }
        free(xhat);
        ewcr_free_wave_csv(&w);
    }
    if (csv) fclose(csv);
    if (man) fclose(man);
    printf("\nRan %d field cases, %d failed.\n", nrun, nfail);
    if (nfail == 0) {
        printf("Field-data five-codec verification PASSED.\n");
        return 0;
    }
    printf("Field-data five-codec verification FAILED.\n");
    return 1;
}
