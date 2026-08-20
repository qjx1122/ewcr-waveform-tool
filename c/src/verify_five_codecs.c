#include "ewcr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int finite_vec(const double *x, int n)
{
    for (int i = 0; i < n; i++) {
        if (!isfinite(x[i])) return 0;
    }
    return 1;
}

static void print_header(void)
{
    printf("\n%-12s %-14s %6s %10s %10s %10s %10s %10s  %s\n",
           "signal", "algo", "N", "CR", "SNR_dB", "PRD", "enc_ms", "dec_ms", "note");
    printf("--------------------------------------------------------------------------------"
           "--------------------------------\n");
}



int main(int argc, char **argv)
{
    int verbose_unit = 1;
    int cycles_short = 2;
    int cycles_long = 4;
    if (argc > 1 && strcmp(argv[1], "--quick") == 0) {
        cycles_long = 2;
    }

    printf("EWCR five-codec C port verification\n");
    printf("===================================\n");
    printf("Unit tests...\n");
    int uf = ewcr_unit_tests(verbose_unit);
    if (uf) {
        fprintf(stderr, "%d unit test(s) failed.\n", uf);
        return 1;
    }
    printf("Unit tests: ALL PASS\n");

    EwcrBenchCfg cfg;
    ewcr_default_cfg(&cfg);

    const char *algos[5] = {"ASBC", "DWT-Hybrid", "CS-OMP", "MMC", "SVDCS"};
    const char *signals_short[4] = {"pure", "sag", "harmonics", "complex"};

    int nmax = (int)round(cfg.fs * cycles_long / cfg.f0) + 16;
    double *x = (double *)ewcr_xmalloc((size_t)nmax * sizeof(double));
    double *xhat = (double *)ewcr_xmalloc((size_t)nmax * sizeof(double));

    mkdir("results", 0755);
    mkdir("c/results", 0755);
    const char *csvpath = "results/verify_five_codecs.csv";
    FILE *csv = fopen(csvpath, "w");
    if (!csv) csv = fopen("c/results/verify_five_codecs.csv", "w");
    if (!csv) csv = fopen("/tmp/verify_five_codecs.csv", "w");
    if (csv) {
        fprintf(csv, "signal,algo,cycles,N,ok,CR,SNR_dB,NMSE_dB,RMSE,PRD,MAXE,PSNR_dB,"
                     "enc_s,dec_s,note\n");
    }

    int nfail = 0, nrun = 0;
    print_header();

    /* --- smoke: 2 cycles x 4 signals x 5 algos --- */
    for (int s = 0; s < 4; s++) {
        int n = 0;
        if (ewcr_gen_ieee1159(signals_short[s], cfg.fs, cycles_short, x, &n) != 0) {
            fprintf(stderr, "signal gen failed: %s\n", signals_short[s]);
            nfail++;
            continue;
        }
        for (int a = 0; a < 5; a++) {
            memset(xhat, 0, (size_t)n * sizeof(double));
            EwcrRunResult r;
            int rc = ewcr_run_codec(algos[a], x, n, &cfg, xhat, &r);
            strncpy(r.signal, signals_short[s], sizeof(r.signal) - 1);
            nrun++;
            int ok = (rc == 0 && r.ok && finite_vec(xhat, r.N > 0 ? r.N : n) &&
                      isfinite(r.metrics.CR) && r.metrics.CR > 0 && isfinite(r.metrics.SNR_dB));
            /* reconstruction quality floor: pure should be reasonably good */
            if (strcmp(signals_short[s], "pure") == 0 && r.metrics.SNR_dB < 12.0) ok = 0;
            if (strcmp(algos[a], "SVDCS") == 0 && strstr(r.note, "LZW roundtrip failed")) ok = 0;
            r.ok = ok;
            if (!ok) nfail++;
            printf("%-12s %-14s %6d %10.3f %10.2f %10.4f %10.2f %10.2f  %s (%s)\n",
                   signals_short[s], algos[a], r.N, r.metrics.CR, r.metrics.SNR_dB, r.metrics.PRD,
                   1000.0 * r.enc_s, 1000.0 * r.dec_s, ok ? "PASS" : "FAIL", r.note);
            if (csv) {
                fprintf(csv, "%s,%s,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,\"%s\"\n",
                        signals_short[s], algos[a], cycles_short, r.N, ok, r.metrics.CR,
                        r.metrics.SNR_dB, r.metrics.NMSE_dB, r.metrics.RMSE, r.metrics.PRD,
                        r.metrics.MAXE, r.metrics.PSNR_dB, r.enc_s, r.dec_s, r.note);
            }
            fflush(stdout);
        }
    }

    /* --- longer pure record (4 cycles) for all five --- */
    {
        int n = 0;
        ewcr_gen_ieee1159("pure", cfg.fs, cycles_long, x, &n);
        printf("\nLonger pure sinusoid (%d cycles, N=%d):\n", cycles_long, n);
        print_header();
        for (int a = 0; a < 5; a++) {
            memset(xhat, 0, (size_t)n * sizeof(double));
            EwcrRunResult r;
            int rc = ewcr_run_codec(algos[a], x, n, &cfg, xhat, &r);
            nrun++;
            int ok = (rc == 0 && r.ok && finite_vec(xhat, r.N > 0 ? r.N : n) &&
                      isfinite(r.metrics.SNR_dB) && r.metrics.SNR_dB >= 12.0 && r.metrics.CR > 0);
            r.ok = ok;
            if (!ok) nfail++;
            printf("%-12s %-14s %6d %10.3f %10.2f %10.4f %10.2f %10.2f  %s (%s)\n", "pure*",
                   algos[a], r.N, r.metrics.CR, r.metrics.SNR_dB, r.metrics.PRD, 1000.0 * r.enc_s,
                   1000.0 * r.dec_s, ok ? "PASS" : "FAIL", r.note);
            if (csv) {
                fprintf(csv, "%s,%s,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,\"%s\"\n",
                        "pure", algos[a], cycles_long, r.N, ok, r.metrics.CR, r.metrics.SNR_dB,
                        r.metrics.NMSE_dB, r.metrics.RMSE, r.metrics.PRD, r.metrics.MAXE,
                        r.metrics.PSNR_dB, r.enc_s, r.dec_s, r.note);
            }
            fflush(stdout);
        }
    }

    if (csv) fclose(csv);
    free(x);
    free(xhat);

    printf("\n===================================\n");
    printf("Ran %d codec cases, %d failed.\n", nrun, nfail);
    if (nfail == 0) {
        printf("Five-codec C verification PASSED.\n");
        return 0;
    }
    printf("Five-codec C verification FAILED.\n");
    return 1;
}
