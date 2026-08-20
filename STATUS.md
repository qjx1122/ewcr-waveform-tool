# STATUS.md
## 当前目标
- 对照 MATLAB 全工程，使 C 五算法功能与 MATLAB 一致（本轮完成对照 + 补齐可选支路）

## 已完成
- [x] MATLAB 工程功能地图：五编码器 + `run_codec_unified` + `unified_metrics` + IEEE 1159 + exp1–exp6
- [x] C 与 MATLAB exp1（10 周期）对照：ASBC/MMC/SVDCS 的 CR/SNR 对齐；CS-OMP CR 对齐
- [x] CS-OMP 增加 MATLAB `basis='dct'`；SVDCS 增加 `sync_mode='zero_cross'`
- [x] `make test` 45/45 PASS

## 进行中
- 无

## 下一步（TODO）
1. 如需：独立 `.ewcr` 解码 CLI
2. DWT `'sym'` 延拓若要对齐 MATLAB CR，需一份能完美重构的 MATLAB `wavedec` 复现
3. 如需：移植 exp2–exp6 扫频（非编解码功能本身）

## 决策记录 / 踩坑
- 功能一致以 `run_codec_unified` + `benchmark_config` 固定工作点为准（MATLAB 实验入口），不是各 `*_encode.m` 的独立默认值（例如 ASBC 文件默认 W=10，基准为 W=20）。
- DWT 周期延拓可完美重构；按 MATLAB `dwt.m`/`idwt.m` 手写 `'sym'` 时 PR 失败，故保持周期延拓，CR 不与 MATLAB 表逐位比较。
- CS-OMP 的 Φ 无法在无 MATLAB ziggurat 表的情况下 bit-identical；算法（高斯测量 + OMP + DFT/DCT 字典）一致，CR 由公式保证一致。
- exp2–exp6 是扫频/资源/作图，不是编解码闭环；C 覆盖五算法 encode→decode。

## 关键文件路径
- `matlab/algorithms/*`、`matlab/common/run_codec_unified.m`、`matlab/common/benchmark_config.m`
- `c/src/{asbc,dwt_hybrid,cs_omp,mmc,svdcs}.c`
- `matlab/results/exp1_fixed_scenarios.csv`
- `c/results/verify_five_codecs.csv`
