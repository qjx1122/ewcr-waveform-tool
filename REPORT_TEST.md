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
