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
