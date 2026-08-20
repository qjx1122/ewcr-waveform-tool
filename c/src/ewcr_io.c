#include "ewcr.h"

#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

void ewcr_buf_init(EwcrBuf *b)
{
    memset(b, 0, sizeof(*b));
}

void ewcr_buf_free(EwcrBuf *b)
{
    free(b->d);
    memset(b, 0, sizeof(*b));
}

static void buf_need(EwcrBuf *b, size_t add)
{
    if (b->n + add <= b->cap) return;
    size_t nc = b->cap ? b->cap : 256;
    while (nc < b->n + add) nc *= 2;
    uint8_t *p = (uint8_t *)realloc(b->d, nc);
    if (!p) {
        fprintf(stderr, "ewcr_buf realloc failed\n");
        exit(2);
    }
    b->d = p;
    b->cap = nc;
}

void ewcr_buf_bytes(EwcrBuf *b, const void *p, size_t n)
{
    buf_need(b, n);
    memcpy(b->d + b->n, p, n);
    b->n += n;
}

void ewcr_buf_u8(EwcrBuf *b, uint8_t v) { ewcr_buf_bytes(b, &v, 1); }

void ewcr_buf_u32(EwcrBuf *b, uint32_t v)
{
    uint8_t t[4] = {(uint8_t)v, (uint8_t)(v >> 8), (uint8_t)(v >> 16), (uint8_t)(v >> 24)};
    ewcr_buf_bytes(b, t, 4);
}

void ewcr_buf_i32(EwcrBuf *b, int32_t v) { ewcr_buf_u32(b, (uint32_t)v); }

void ewcr_buf_f64(EwcrBuf *b, double v)
{
    uint64_t u;
    memcpy(&u, &v, 8);
    uint8_t t[8];
    for (int i = 0; i < 8; i++) t[i] = (uint8_t)((u >> (8 * i)) & 0xFF);
    ewcr_buf_bytes(b, t, 8);
}

void ewcr_buf_str(EwcrBuf *b, const char *s)
{
    uint32_t n = s ? (uint32_t)strlen(s) : 0;
    ewcr_buf_u32(b, n);
    if (n) ewcr_buf_bytes(b, s, n);
}

void ewcr_pack_header(EwcrBuf *b, uint8_t algo_id, int n, double fs, double f0, double bits_compressed)
{
    ewcr_buf_bytes(b, "EWCR", 4);
    ewcr_buf_u8(b, 1); /* version */
    ewcr_buf_u8(b, algo_id);
    ewcr_buf_u8(b, 0);
    ewcr_buf_u8(b, 0);
    ewcr_buf_u32(b, (uint32_t)n);
    ewcr_buf_f64(b, fs);
    ewcr_buf_f64(b, f0);
    ewcr_buf_f64(b, bits_compressed);
}

uint8_t ewcr_algo_id(const char *algo)
{
    if (!algo) return 0;
    if (strcmp(algo, "ASBC") == 0) return 1;
    if (strcmp(algo, "DWT-Hybrid") == 0 || strcmp(algo, "DWT") == 0) return 2;
    if (strcmp(algo, "CS-OMP") == 0 || strcmp(algo, "CS") == 0) return 3;
    if (strcmp(algo, "MMC") == 0) return 4;
    if (strcmp(algo, "SVDCS") == 0) return 5;
    return 0;
}

int ewcr_write_wave_csv(const char *path, const double *x, int n, double fs)
{
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;
    fprintf(fp, "sample,time_s,value\n");
    for (int i = 0; i < n; i++) {
        double t = (fs > 0) ? ((double)i / fs) : (double)i;
        fprintf(fp, "%d,%.9g,%.9g\n", i, t, x[i]);
    }
    fclose(fp);
    return 0;
}

int ewcr_write_compare_csv(const char *path, const double *x, const double *xhat, int n, double fs)
{
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;
    fprintf(fp, "sample,time_s,original,reconstructed,error\n");
    for (int i = 0; i < n; i++) {
        double t = (fs > 0) ? ((double)i / fs) : (double)i;
        fprintf(fp, "%d,%.9g,%.9g,%.9g,%.9g\n", i, t, x[i], xhat[i], x[i] - xhat[i]);
    }
    fclose(fp);
    return 0;
}

int ewcr_write_bytes(const char *path, const uint8_t *data, int nbytes)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;
    if (nbytes > 0 && data) fwrite(data, 1, (size_t)nbytes, fp);
    fclose(fp);
    return 0;
}

int ewcr_save_case_files(const char *dir, const char *stem, const double *x, const double *xhat,
                         int n, double fs, const EwcrRunResult *r)
{
    if (!dir || !stem || n <= 0) return -1;
    mkdir(dir, 0755);
    char p[768];
    snprintf(p, sizeof(p), "%s/%s_original.csv", dir, stem);
    if (ewcr_write_wave_csv(p, x, n, fs) != 0) return -2;
    snprintf(p, sizeof(p), "%s/%s_reconstructed.csv", dir, stem);
    if (ewcr_write_wave_csv(p, xhat, n, fs) != 0) return -3;
    snprintf(p, sizeof(p), "%s/%s_compare.csv", dir, stem);
    ewcr_write_compare_csv(p, x, xhat, n, fs);
    if (r && r->compressed && r->compressed_nbytes > 0) {
        snprintf(p, sizeof(p), "%s/%s_compressed.ewcr", dir, stem);
        if (ewcr_write_bytes(p, r->compressed, r->compressed_nbytes) != 0) return -4;
    }
    return 0;
}

void ewcr_result_release(EwcrRunResult *r)
{
    if (!r) return;
    free(r->compressed);
    r->compressed = NULL;
    r->compressed_nbytes = 0;
}
