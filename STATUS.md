# STATUS.md
## 当前目标
- 对照 `matlab/` 复查 C 五算法功能是否正常（已完成）

## 已完成
- [x] 逐文件对照 MATLAB 五编码器 / 默认配置 / `unified_metrics` / IEEE 1159 与 C 实现
- [x] 修复 ASBC `interpft` 偶数长度 Nyquist 分裂（MATLAB `ceil((m+1)/2)`）
- [x] 10 周期 exp1 对照：ASBC/MMC/SVDCS 的 CR 与 MATLAB 表一致；MMC/SVDCS SNR 一致
- [x] `make test` 45/45 PASS；`make test-data` 120/120 PASS

## 进行中
- 无

## 下一步（TODO）
1. 如需：独立 `.ewcr` 解码 CLI
2. 如需：DWT 改为 MATLAB `sym` 延拓以对齐 CR
3. 如需：整段长记录（>10 周期）

## 决策记录 / 踩坑
- `ewcr_interpft` 曾用整数 `(n+1)/2`（floor），偶数长度时把 DC 当 Nyquist 对半切开；ASBC 4 周期纯正弦 SNR 从 ~304 dB 掉到 ~102 dB。改为 `(n+2)/2` = MATLAB `ceil((n+1)/2)` 后恢复。
- 10 周期纯正弦 ASBC SNR C=266 dB / MATLAB=294 dB，均为数值无损；CR 精确一致。差来自 1280 点 Bluestein vs MATLAB FFT 舍入。
- DWT 周期延拓 vs MATLAB `sym`：CR 不可比，C 闭环 SNR 反而更高。
- CS-OMP Φ RNG 不同，CR 公式一致故 CR 必等；SNR 只应接近。
- 现场 MMC 窗 N=200 非 2 幂，C 跳过 Haar 残差（MATLAB Haar 在非 2 幂上也不合法）。

## 关键文件路径
- `matlab/algorithms/*` 与 `c/src/{asbc,dwt_hybrid,cs_omp,mmc,svdcs}.c`
- `matlab/results/exp1_fixed_scenarios.csv`
- `c/results/verify_five_codecs.csv`
- `c/results/verify_field_data.csv`
