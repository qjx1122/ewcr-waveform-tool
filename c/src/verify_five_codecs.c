#include "ewcr.h"

#include <math.h>
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
           "signal", "algo", "N", "压缩比CR", "相似度%", "SNR_dB", "enc_ms", "dec_ms", "note");
    printf("--------------------------------------------------------------------------------"
           "--------------------------------\n");
}



int main(int argc, char **argv)
{
    int verbose_unit = 1;
    int cycles_short = 2;
    int cycles_long = 4;
    int do_matlab = 1;
    if (argc > 1 && strcmp(argv[1], "--quick") == 0) {
        cycles_long = 2;
        do_matlab = 0;
    }

    printf("EWCR five-codec C port verification\n");
    printf("===================================\n");
    printf("压缩比 CR = (N*16 bit)/bits_compressed;  相似度%% = 100*Pearson(x,xhat)\n");
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

    const int cycles_matlab = 10; /* matlab/experiments/exp1_fixed_scenarios.m */
    int nmax = (int)round(cfg.fs * cycles_matlab / cfg.f0) + 16;
    double *x = (double *)ewcr_xmalloc((size_t)nmax * sizeof(double));
    double *xhat = (double *)ewcr_xmalloc((size_t)nmax * sizeof(double));

    mkdir("results", 0755);
    mkdir("results/out", 0755);
    mkdir("c/results", 0755);
    const char *csvpath = "results/verify_five_codecs.csv";
    FILE *csv = fopen(csvpath, "w");
    if (!csv) csv = fopen("c/results/verify_five_codecs.csv", "w");
    if (!csv) csv = fopen("/tmp/verify_five_codecs.csv", "w");
    if (csv) {
        fprintf(csv, "signal,algo,cycles,N,ok,compress_ratio,similarity_pct,corr,SNR_dB,"
                     "NMSE_dB,RMSE,PRD,MAXE,PSNR_dB,enc_s,dec_s,note\n");
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
            printf("%-12s %-14s %6d %10.3f %10.4f %10.2f %10.2f %10.2f  %s (%s)\n",
                   signals_short[s], algos[a], r.N, r.metrics.CR, r.metrics.similarity_pct,
                   r.metrics.SNR_dB, 1000.0 * r.enc_s, 1000.0 * r.dec_s, ok ? "PASS" : "FAIL",
                   r.note);
            if (csv) {
                fprintf(csv, "%s,%s,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,\"%s\"\n",
                        signals_short[s], algos[a], cycles_short, r.N, ok, r.metrics.CR,
                        r.metrics.similarity_pct, r.metrics.corr, r.metrics.SNR_dB,
                        r.metrics.NMSE_dB, r.metrics.RMSE, r.metrics.PRD, r.metrics.MAXE,
                        r.metrics.PSNR_dB, r.enc_s, r.dec_s, r.note);
            }
            if (ok) {
                char algo_id[32];
                snprintf(algo_id, sizeof(algo_id), "%s", algos[a]);
                for (char *p = algo_id; *p; p++) if (*p == '-') *p = '_';
                char stem[128];
                snprintf(stem, sizeof(stem), "synth_%s_%s", signals_short[s], algo_id);
                int nn = r.N > 0 ? r.N : n;
                if (nn > n) nn = n;
                ewcr_save_case_files("results/out", stem, x, xhat, nn, cfg.fs, &r);
            }
            ewcr_result_release(&r);
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
            printf("%-12s %-14s %6d %10.3f %10.4f %10.2f %10.2f %10.2f  %s (%s)\n", "pure*",
                   algos[a], r.N, r.metrics.CR, r.metrics.similarity_pct, r.metrics.SNR_dB,
                   1000.0 * r.enc_s, 1000.0 * r.dec_s, ok ? "PASS" : "FAIL", r.note);
            if (csv) {
                fprintf(csv, "%s,%s,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,\"%s\"\n",
                        "pure", algos[a], cycles_long, r.N, ok, r.metrics.CR,
                        r.metrics.similarity_pct, r.metrics.corr, r.metrics.SNR_dB,
                        r.metrics.NMSE_dB, r.metrics.RMSE, r.metrics.PRD, r.metrics.MAXE,
                        r.metrics.PSNR_dB, r.enc_s, r.dec_s, r.note);
            }
            fflush(stdout);
        }
    }

    /* --- 10-cycle IEEE 1159 vs MATLAB exp1_fixed_scenarios.csv --- */
    if (do_matlab) {
        typedef struct {
            const char *signal;
            const char *algo;
            double CR;
            double SNR_dB;
        } MatlabRef;
        static const MatlabRef refs[] = {
            {"pure", "ASBC", 86.7796610169491, 294.594527938571},
            {"pure", "DWT-Hybrid", 5.40226853073068, 35.5233084004896},
            {"pure", "CS-OMP", 4.90421455938697, 94.6518489994766},
            {"pure", "MMC", 29.2571428571429, 79.9210537311647},
            {"pure", "SVDCS", 3.19201995012469, 90.3087336228232},
            {"sag", "ASBC", 44.7161572052402, 33.9734530777162},
            {"sag", "DWT-Hybrid", 4.0, 34.6729029763487},
            {"sag", "CS-OMP", 4.90421455938697, 43.2314356453789},
            {"sag", "MMC", 20.3073872087258, 76.5696998844263},
            {"sag", "SVDCS", 3.12004875076173, 72.9193581731085},
            {"harmonics", "ASBC", 60.5917159763314, 33.6968671628163},
            {"harmonics", "DWT-Hybrid", 2.81821934773634, 37.5724165203979},
            {"harmonics", "CS-OMP", 4.90421455938697, 95.2428189778234},
            {"harmonics", "MMC", 16.4234161988773, 52.9568680966022},
            {"harmonics", "SVDCS", 3.02779420461266, 79.3615271150213},
            {"complex", "ASBC", 8.92763731473409, 28.4978610851198},
            {"complex", "DWT-Hybrid", 3.33251972988365, 35.8367707115774},
            {"complex", "CS-OMP", 4.90421455938697, 26.8175234964565},
            {"complex", "MMC", 16.3122262046993, 34.9052431400195},
            {"complex", "SVDCS", 1.82726623840114, 29.0012529042353},
        };
        const char *kinds[4] = {"pure", "sag", "harmonics", "complex"};
        printf("\nMATLAB exp1 compare (10 cycles, N=2560, fs=12800):\n");
        printf("%-12s %-14s %10s %10s %10s %10s %10s %8s\n", "signal", "algo", "CR_C",
               "CR_ML", "SNR_C", "SNR_ML", "sim%", "align");
        printf("--------------------------------------------------------------------------------"
               "----------------\n");
        for (int s = 0; s < 4; s++) {
            int n = 0;
            if (ewcr_gen_ieee1159(kinds[s], cfg.fs, cycles_matlab, x, &n) != 0) {
                fprintf(stderr, "signal gen failed: %s\n", kinds[s]);
                nfail++;
                continue;
            }
            for (int a = 0; a < 5; a++) {
                memset(xhat, 0, (size_t)n * sizeof(double));
                EwcrRunResult r;
                int rc = ewcr_run_codec(algos[a], x, n, &cfg, xhat, &r);
                nrun++;
                const MatlabRef *ref = NULL;
                for (size_t k = 0; k < sizeof(refs) / sizeof(refs[0]); k++) {
                    if (strcmp(refs[k].signal, kinds[s]) == 0 &&
                        strcmp(refs[k].algo, algos[a]) == 0) {
                        ref = &refs[k];
                        break;
                    }
                }
                int ok = (rc == 0 && r.ok && finite_vec(xhat, r.N > 0 ? r.N : n) &&
                          isfinite(r.metrics.CR) && r.metrics.CR > 0 && isfinite(r.metrics.SNR_dB));
                if (strcmp(kinds[s], "pure") == 0 && r.metrics.SNR_dB < 12.0) ok = 0;
                int align = 1;
                const char *why = "ok";
                if (!ref) {
                    align = 0;
                    why = "no-ref";
                } else if (strcmp(algos[a], "CS-OMP") == 0) {
                    /* CR is analytic in (N, M_ratio, ybits); Phi RNG does not affect CR. */
                    if (fabs(r.metrics.CR - ref->CR) / ref->CR > 1e-9) {
                        align = 0;
                        why = "CR";
                    }
                } else if (strcmp(algos[a], "MMC") == 0) {
                    if (fabs(r.metrics.CR - ref->CR) / ref->CR > 1e-6) {
                        align = 0;
                        why = "CR";
                    }
                    if (strcmp(kinds[s], "pure") == 0 &&
                        fabs(r.metrics.SNR_dB - ref->SNR_dB) > 0.05) {
                        align = 0;
                        why = "SNR";
                    }
                } else if (strcmp(algos[a], "ASBC") == 0 && strcmp(kinds[s], "pure") == 0) {
                    if (r.metrics.SNR_dB < 200.0) {
                        align = 0;
                        why = "SNR";
                    }
                } else if (strcmp(algos[a], "SVDCS") == 0 && strcmp(kinds[s], "pure") == 0) {
                    if (fabs(r.metrics.SNR_dB - ref->SNR_dB) > 1.0) {
                        align = 0;
                        why = "SNR";
                    }
                    if (fabs(r.metrics.CR - ref->CR) / ref->CR > 0.05) {
                        align = 0;
                        why = "CR";
                    }
                } else if (strcmp(algos[a], "DWT-Hybrid") == 0) {
                    why = "ext";
                    if (strcmp(kinds[s], "pure") == 0 && r.metrics.SNR_dB < 20.0) {
                        align = 0;
                        why = "SNR";
                    }
                }
                if (!ok || !align) {
                    nfail++;
                    r.ok = 0;
                }
                printf("%-12s %-14s %10.4f %10.4f %10.2f %10.2f %10.4f %8s\n", kinds[s], algos[a],
                       r.metrics.CR, ref ? ref->CR : NAN, r.metrics.SNR_dB,
                       ref ? ref->SNR_dB : NAN, r.metrics.similarity_pct,
                       (!ok || !align) ? "FAIL" : why);
                if (csv) {
                    fprintf(csv,
                            "%s,%s,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,"
                            "\"matlab-exp1 %s\"\n",
                            kinds[s], algos[a], cycles_matlab, r.N, ok && align, r.metrics.CR,
                            r.metrics.similarity_pct, r.metrics.corr, r.metrics.SNR_dB,
                            r.metrics.NMSE_dB, r.metrics.RMSE, r.metrics.PRD, r.metrics.MAXE,
                            r.metrics.PSNR_dB, r.enc_s, r.dec_s, why);
                }
                ewcr_result_release(&r);
                fflush(stdout);
            }
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
