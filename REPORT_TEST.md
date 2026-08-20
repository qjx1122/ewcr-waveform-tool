# REPORT_TEST.md

## [2026-08-20] 专题：MATLAB 五算法 → C 移植与验证
- 类型：验证专题
- 目标与假设：
  - 将 `matlab/algorithms` 中 ASBC / DWT-Hybrid / CS-OMP / MMC / SVDCS 编解码链路移植为可独立编译的 C11 实现
  - 用 IEEE 1159 合成波形做 encode→decode 闭环，确认有限重构、正压缩比、纯正弦 SNR≥12 dB
- 方法 / 数据 / 参数：
  - 代码：`c/src/*.c`，入口 `c/src/verify_five_codecs.c`
  - 信号：`pure` / `sag` / `harmonics` / `complex`，fs=12800 Hz，2 周期（N=512）+ 纯正弦 4 周期（N=1024）
  - 参数对齐 `matlab/common/benchmark_config.m` 固定工作点
  - 内核单测：FFT、DCT、Haar、db4 PR、Huffman、最小二乘
- 结果 / 结论：
  - 8/8 单元测试通过；25/25 编解码用例通过（`make test`）
  - 纯正弦 2 周期：ASBC SNR≈304 dB；CS-OMP ≈91 dB；MMC CR=29.257、SNR=79.92 dB（与 MATLAB `exp1` 纯正弦参考一致）；SVDCS ≈90 dB；DWT-Hybrid ≈45 dB
  - 短记录上 DWT/SVDCS 的 CR 受 Huffman / MC+LZW 头开销影响，可 <1；属记账与帧数问题，不是解码错误
  - CS-OMP 的 Φ 使用本地 xorshift+Box-Muller，与 MATLAB `rng(42)+randn` 不同，SNR 不可逐点对比
- 是否进入 REPORT.md（稳定结论）：否
- 遗留问题：
  - 未把 exp1–exp6 全套扫频搬到 C
  - DWT 使用周期延拓 db4，不是 MATLAB `wavedec` 默认 `sym` 延拓
  - 本环境无 MATLAB，无法做 C vs MATLAB 逐样本 bit-identical 对照

## [2026-08-20] 专题：现场 CSV（data/）五算法验证
- 类型：验证专题
- 目标与假设：用仓库 `data/wave1..4.csv` 真实电气波形，对 C 移植的五算法做 encode→decode 闭环
- 方法 / 数据 / 参数：
  - 数据：10 kHz，约 5–10 s，列 timestamp,UA,IA,UB,IB,UC,IC；电压峰值约 330 V，电流约 1 A
  - 片段：跳过第 1 个 50 Hz 周期后取 10 周期（N=2000）的 UA 与 IA
  - 入口：`c/src/verify_field_data.c`，`make test-data`
  - MMC 窗 N=200（非 2 幂）只用 none/DCT 残差；SVDCS `signal_scale`=片段峰值
- 结果 / 结论：
  - 40/40 通过（4 文件 × 2 通道 × 5 算法）
  - 电压 UA：MMC CR≈16.6–16.9、SNR≈43–46 dB；DWT SNR≈48–49 dB；ASBC CR≈11–15、SNR≈27 dB；CS-OMP CR=4.88、SNR≈43–46 dB；SVDCS SNR≈39–46 dB
  - 电流 IA：SNR 普遍低于电压（ASBC 约 15–24 dB，其余约 28–45 dB），CR 模式与电压类似
- 是否进入 REPORT.md（稳定结论）：否
- 遗留问题：未覆盖 B/C 相与整段长记录；ASBC F0 栅格偏差

## [2026-08-20] 专题：压缩比与波形相似度输出
- 类型：验证专题
- 目标与假设：在原波 vs 复原波对比中显式输出压缩比与相似度
- 方法 / 数据 / 参数：
  - `compress_ratio = CR = (N×16)/bits_compressed`
  - `corr = Pearson(x, xhat)`，`similarity_pct = 100×corr`
- 结果 / 结论：
  - 现场 UA：压缩比 MMC≈16.6、CS-OMP=4.88、ASBC≈11–15；相似度均 >99.9%
  - Pearson 对工频波形极不敏感（小误差仍 ~100%），细部差异仍看 SNR_dB
- 是否进入 REPORT.md（稳定结论）：否
- 遗留问题：未输出逐点对比波形文件（已在后续专题解决）

## [2026-08-20] 专题：压缩文件与复原文件落盘
- 类型：用户专题
- 目标与假设：每次编解码写出压缩码流和复原波形，便于离线查看
- 方法 / 数据 / 参数：`ewcr_save_case_files`；现场 10 周期 UA/IA；合成 IEEE 1159
- 结果 / 结论：
  - `c/results/out/` 每例 4 个文件：original / reconstructed / compare CSV + compressed.ewcr
  - 例：`wave1_UA_MMC_compressed.ewcr` 约 4.5 KiB，对应 2000 点电压窗
- 是否进入 REPORT.md（稳定结论）：否
- 遗留问题：无独立 `ewcr_decode` 命令行

## [2026-08-20] 专题：data/ 全通道重跑
- 类型：验证专题
- 目标与假设：wave1–4 的 UA/IA/UB/IB/UC/IC 全部进入五算法验证
- 方法 / 数据 / 参数：10 周期 @ 10 kHz，跳过首周期
- 结果 / 结论：120/120 PASS。电压通道相似度普遍 >99.9%；wave4 IB 偏置电流相似度可低至 ~63%（CS-OMP）
- 是否进入 REPORT.md（稳定结论）：否
- 遗留问题：未覆盖整段长记录

## [2026-08-20] 专题：对照 MATLAB 复查 C 功能
- 类型：验证专题
- 目标与假设：C 移植应在相同工作点上复现 MATLAB 五算法的压缩比与重构保真度，而不仅是 C 编码器/解码器自洽
- 方法 / 数据 / 参数：
  - 对照 `matlab/algorithms/{ASBC,DWT_Hybrid,CS_OMP,MMC,SVDCS}` 与 `c/src/*.c`
  - 金标准：`matlab/results/exp1_fixed_scenarios.csv`（IEEE 1159，10 周期，fs=12800，`benchmark_config` 固定点）
  - 入口：`c/src/verify_five_codecs.c` 的 MATLAB exp1 compare 段；内核单测含 naive DFT 与 MATLAB `interpft`
- 结果 / 结论：
  - **缺陷**：`ewcr_interpft` 用整数 `(n+1)/2` 代替 MATLAB `ceil((n+1)/2)`，偶数长度把 DC 当 Nyquist 对半切。ASBC 4 周期纯正弦 SNR 从 ~304 dB 掉到 ~102 dB。修复后 4 周期 SNR≈304 dB，10 周期 ASBC CR 与 MATLAB **精确一致**（pure 86.7797，sag 44.7162，harmonics 60.5917，complex 8.9276）
  - **MMC**：四场景 CR/SNR 与 MATLAB 表一致（pure CR=29.2571 SNR=79.92 dB）
  - **SVDCS**：四场景 CR/SNR 与 MATLAB 表一致（pure CR=3.1920 SNR=90.31 dB）
  - **CS-OMP**：CR 四场景均为 4.9042（公式量）；SNR 接近但因 Φ 的 RNG 不同不完全相等
  - **DWT-Hybrid**：周期延拓 db4 ≠ MATLAB `wavedec` 默认 `sym`，CR 不可比；C 闭环 SNR 约 39–45 dB，功能正常
  - 现场 `data/` 120/120 仍 PASS
- 是否进入 REPORT.md（稳定结论）：否
- 遗留问题：纯正弦 ASBC 10 周期 SNR 266 vs 294 dB（均为数值无损，CR 已对齐）；DWT 若要对齐 MATLAB CR 需改延拓

## [2026-08-20] 专题：MATLAB 全工程功能 vs C 对齐
- 类型：验证专题
- 目标与假设：C 应实现 MATLAB 五算法编解码功能（同一工作点、同一码率公式、同一重构链路）
- 方法 / 数据 / 参数：
  - MATLAB 入口 `run_all_experiments` → exp1 固定场景为金标准
  - C：`make test` 10 周期对照 `exp1_fixed_scenarios.csv`
- 结果 / 结论：
  - MATLAB 功能 = 五编码器 + 统一指标 + IEEE 1159 + 六组实验（RD/时长/噪声/耗时/载荷）。C 覆盖五编码器闭环与 exp1 工作点，不覆盖 exp2–exp6 扫频作图。
  - ASBC/MMC/SVDCS：CR 与 MATLAB 表精确一致；MMC/SVDCS SNR 精确一致；ASBC 扰动场景 SNR 一致，纯正弦均为数值无损。
  - CS-OMP：CR 精确一致；已实现 DFT（默认）与 DCT 字典；Φ 的 RNG 不同故 SNR 接近而非相同。
  - DWT-Hybrid：阶段与 MATLAB 相同（db4→阈值→量化→delta→Huffman）；延拓为周期而非 `'sym'`，CR 不同，C 闭环 SNR 可用。
  - SVDCS 已实现 `'nominal'` 与 `'zero_cross'` 同步。
- 是否进入 REPORT.md（稳定结论）：否
- 遗留问题：无 MATLAB Wavelet Toolbox 环境下的 `'sym'` 逐样本金标；CS-OMP 无法复现 MATLAB ziggurat `randn`
