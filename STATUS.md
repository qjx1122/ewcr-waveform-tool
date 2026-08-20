# STATUS.md
## 当前目标
- 将 `matlab/` 五算法编解码移植为 C 并完成验证测试（已完成本会话范围）

## 已完成
- [x] 按 BOOTSTRAP 建立会话上下文
- [x] C11 移植五算法：ASBC、DWT-Hybrid、CS-OMP、MMC、SVDCS
- [x] IEEE 1159 合成信号、FFT/DCT/Haar/db4、Huffman、12-bit LZW
- [x] `c/build/verify_five_codecs`：8 项内核单测 + 25 项编解码用例全部通过
- [x] 结果表 `c/results/verify_five_codecs.csv`

## 进行中
- 无

## 下一步（TODO）
1. 如需：把 exp1–exp6 扫频/噪声/资源实验也做到 C
2. 如有 MATLAB：做 C vs MATLAB 逐样本对照（CS 的 Φ 需对齐 RNG）
3. 短记录上 DWT/SVDCS 的 CR 头开销可再优化记账或打包

## 决策记录 / 踩坑
- Arena 分支固定 `arena/01a01e0b-ewcr-waveform-tool`，不另开 feature 分支
- 沙箱无 MATLAB/FFTW 头文件：自实现混合基 FFT（pow2 + Bluestein）
- 周期延拓 db4 必须用同一套正交滤波器、卷积索引 `(2i-k) mod n` 才能 PR
- CS-OMP 的 A=ΦΨ 用 Φ 各行 FFT 构造，避免显式 N×N DFT 矩阵
- MMC 纯正弦 CR=29.257、SNR=79.92 dB，与 MATLAB `exp1_fixed_scenarios.csv` 纯正弦行一致
- 验证口径：有限重构、CR>0、纯正弦 SNR≥12 dB；SVDCS 额外检查 LZW 往返

## 关键文件路径
- `c/README.md` — C 工程说明
- `c/Makefile` — `make` / `make test`
- `c/include/ewcr.h`
- `c/src/{asbc,dwt_hybrid,cs_omp,mmc,svdcs,ewcr_math,ewcr_signals,verify_five_codecs}.c`
- `c/results/verify_five_codecs.csv`
- `matlab/` — 原始 MATLAB 基准
- `BOOTSTRAP.md` / `STATUS.md` / `session/NILM_AC_session_complete.md` / `REPORT_TEST.md`
