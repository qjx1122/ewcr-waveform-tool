# Five-codec electrical waveform compression — C port

Self-contained C11 translation of the MATLAB codecs under `matlab/algorithms/`:

| Codec | MATLAB | C |
| --- | --- | --- |
| ASBC | `asbc_encode/decode` | `src/asbc.c` |
| DWT-Hybrid | `dwt_hybrid_encode/decode` | `src/dwt_hybrid.c` |
| CS-OMP | `cs_omp_encode/decode` | `src/cs_omp.c` |
| MMC | `mmc_encode_signal/decode_signal` | `src/mmc.c` |
| SVDCS | `svdcs_encode/decode` + LZW | `src/svdcs.c` |

IEEE Std 1159 test signals, FFT/DCT/Haar/db4 DWT, Huffman and 12-bit LZW live in `src/ewcr_math.c` and `src/ewcr_signals.c`. No MATLAB, FFTW or BLAS required (`-lm` only).

## Build and verify

```bash
cd c
make
make test          # IEEE 1159 synthetic: ./build/verify_five_codecs
make test-data     # field CSV in ../data: ./build/verify_field_data
```

`make test` runs:

1. Kernel unit tests (FFT vs naive DFT, MATLAB `interpft` even/odd, DCT, Haar, db4 PR, Huffman, least squares)
2. Encode→decode for all **five algorithms** on IEEE 1159 `pure` / `sag` / `harmonics` / `complex` (2 cycles @ 12.8 kHz)
3. A longer pure sinusoid (4 cycles) for all five codecs
4. 10-cycle compare against `matlab/results/exp1_fixed_scenarios.csv` (same `benchmark_config` operating point)

`make test-data` scans every `*.csv` in `data/` (currently wave1–4, 10 kHz three-phase), takes 10 cycles of **all six channels UA/IA/UB/IB/UC/IC** after skipping the first cycle, and runs the same five codecs (4 files × 6 channels × 5 algos = 120 cases). SVDCS `signal_scale` is set to the segment peak so 16-bit p.u. quantizers do not clip field voltages.

CSV: `c/results/verify_five_codecs.csv`, `c/results/verify_field_data.csv`

Per-case files (under `c/results/out/`):

- `*_original.csv` — original waveform (`sample,time_s,value`)
- `*_reconstructed.csv` — decoded waveform
- `*_compare.csv` — original vs reconstructed vs error
- `*_compressed.ewcr` — binary compressed bitstream (magic `EWCR`)
- `manifest.csv` — index of field-data outputs

A case PASSes when reconstruction is finite, `CR > 0`, `SNR` is finite, and SNR on a pure sinusoid is at least 12 dB. SVDCS also checks the LZW payload round-trip.

Primary columns (original `x` vs reconstructed `xhat`):

- **压缩比 `compress_ratio` / `CR`** = `(N × 16 bit) / bits_compressed`
- **波形相似度 `similarity_pct`** = `100 × PearsonCorr(x, xhat)`（%）；`corr` 为同一系数，范围 [-1, 1]
- `SNR_dB` remains the energy-fidelity metric `10 log10(||x||² / ||x-xhat||²)`

## Defaults (aligned with `matlab/common/benchmark_config.m`)

- `fs = 12800`, `f0 = 50`, 16-bit original PCM
- ASBC: `W=20`, `Rk=8`, 5-cycle blocks
- DWT-Hybrid: periodized db4, energy retain 0.99, 8-bit quant, delta + Huffman
- CS-OMP: `M/N=0.2`, unitary DFT dictionary, `K0=48` (Phi is a local RNG, not MATLAB `mt19937ar`)
- MMC: window `N = fs/f0`, `n_tot = N` (1 bit/sample), same model/residual search as MATLAB
- SVDCS: nominal-cycle Lagrange sync, `G=0.03`, `β=0.1`, 16-bit FFT bins, 12-bit LZW

## MATLAB vs C

10-cycle IEEE 1159 @ 12.8 kHz (`N=2560`), C vs `matlab/results/exp1_fixed_scenarios.csv`:

| signal | algo | CR C | CR MATLAB | SNR C | SNR MATLAB |
| --- | --- | ---: | ---: | ---: | ---: |
| pure | ASBC | 86.7797 | 86.7797 | 266.20 | 294.59 |
| pure | MMC | 29.2571 | 29.2571 | 79.92 | 79.92 |
| pure | SVDCS | 3.1920 | 3.1920 | 90.31 | 90.31 |
| pure | CS-OMP | 4.9042 | 4.9042 | 94.98 | 94.65 |
| sag / harmonics / complex | ASBC, MMC, SVDCS | match to printed digits | | match to printed digits | |

- **ASBC / MMC / SVDCS** share the MATLAB architecture, bit accounting, and (for MMC/SVDCS) reconstruction SNR on these records.
- **CS-OMP** CR is analytic in `(N, M/N, ybits)` and matches; Φ uses a local xorshift+Box-Muller, not MATLAB `rng(42)+randn`, so SNR is close but not identical.
- **DWT-Hybrid** follows the same stages (db4 DWT → energy/universal threshold → quant → delta → Huffman) with **periodized** db4 (perfect reconstruction). MATLAB `wavedec` default `dwtmode` is `'sym'`, so DWT CR is not identical; C encode/decode is self-consistent.
- CS-OMP also implements MATLAB `basis='dct'` (`dctmtx`); the benchmark default remains `'dft'`. Φ still uses a local RNG, not MATLAB `mt19937ar` + ziggurat `randn`.
- SVDCS implements MATLAB `sync_mode` `'nominal'` (benchmark) and `'zero_cross'`.
- `interpft` (ASBC baseband upsample) uses MATLAB `ceil((m+1)/2)` Nyquist splitting.
