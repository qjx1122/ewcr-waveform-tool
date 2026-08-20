# Five-Codes Electrical Waveform Compression Benchmark - MATLAB

本工程统一比较 5 种电气波形压缩/重构算法：

1. **ASBC** - Adaptive Subband Compression
2. **DWT-Hybrid** - DWT + adaptive threshold + quantization + delta + Huffman
3. **CS-OMP** - compressed sensing + orthogonal matching pursuit
4. **MMC** - multiple-model coding，受总码率约束的两阶段工程复现
5. **SVDCS** - spectral-variation novelty compression + LZW

v2 在上一版的压缩比、重构误差和 Rate-Distortion 之外，新增了用户关心的三类工程指标：

- **压缩后究竟保留什么参数/码流字段**；
- **编码和解码耗时**：冷启动、稳态中位数、P90、吞吐率、实时因子；
- **边端资源**：逻辑压缩载荷、MATLAB struct/cell 内存、当前 MATLAB 实现的主要工作内存估计、持续状态内存估计、MATLAB 源码体积。

## 推荐环境

- MATLAB R2025b
- Signal Processing Toolbox
- Wavelet Toolbox
- Communications Toolbox

## 一键运行

```matlab
cd('FiveCodesBenchmark_MATLAB');
run_all_experiments;
```

## 新增实验

### Experiment 5 - Runtime / RAM / code footprint

输出：

- `results/exp5_runtime_memory.csv`
- `results/exp5_runtime_memory_scaling.csv`
- `results/exp5_runtime_memory.mat`
- `figures/exp5_runtime.png`
- `figures/exp5_working_memory.png`
- `figures/exp5_payload_memory.png`

时间测试流程：

1. 记录一次 cold run；
2. 做 `cfg.profile.warmup` 次预热并丢弃；
3. 做 `cfg.profile.repeats` 次正式计时；
4. 报告 encoder/decoder 的 median、mean、P90；
5. 同时计算处理时间 / 信号时长的 **real-time factor**。`<1` 表示在该 MATLAB 主机上处理速度快于数据实时到达速度。

### Experiment 6 - Compressed payload manifest

输出：

- `results/exp6_payload_manifest.csv`

该表明确记录每个算法压缩后需要保存/传输的内容，以及哪些字段只是当前 MATLAB 工程为了调试或独立验证而保留、实际部署可去掉。

## 资源字段如何理解

### `logical_payload_bytes`

按算法自身 bit accounting 得到的真实逻辑压缩载荷大小，等于 `bits_compressed/8`。这是比较通信/存储最重要的量。

### `matlab_packet_bytes`

`whos` 统计当前 MATLAB struct/cell 对象在内存中的大小。它通常比实际压缩码流大很多，因为 MATLAB 有容器开销，而且当前工程保留了部分调试/重构字段。

### `encoder_working_est_bytes` / `decoder_working_est_bytes`

根据当前 MATLAB 代码中主要同时驻留的数值矩阵/缓存给出的确定性工程估计。它用于比较算法级 RAM 压力，但**不包含 MATLAB JIT、Toolbox、OS、内存碎片等运行时开销**，也不等于最终 MCU/FPGA RAM。

### `edge_state_est_bytes`

算法跨块/跨周期需要持续保留的状态和基本缓冲区估计。例如 MMC 需要前两个重构窗口；SVDCS 需要 MAX/MIN/SUM/CNT；CS 的当前实现工作矩阵很大，但可由 seed/基函数重新生成，因此持久状态与工作 RAM 是两个概念。

### `source_mfile_bytes`

算法目录下 `.m` 文件总大小，只用于比较 MATLAB 源码体积。它**绝不能当作嵌入式固件 flash 占用**。真实 flash/RAM 必须用 MATLAB Coder/Embedded Coder 或目标编译器生成代码后读 map/report。

## 一个特别重要的公平性说明

- 当前 **CS-OMP** 使用显式 dense `Phi`、`Psi`、`A` 矩阵；DFT 基下 `Psi/A` 为复数 double，因此 RAM 很高。结构化 sensing、FFT 隐式基和边生成矩阵可以把资源需求降低几个数量级，但那属于另一个工程实现，不能与当前 MATLAB 代码混为一谈。
- 当前 **MMC** 编码器是模型候选 + 比特分配 + residual method 的搜索器，因此 encode 比 decode 重得多；解码只需要执行最终选中的一个模型和一个残差解码器。
- 当前 **SVDCS** MATLAB LZW 为便于可逆验证，使用 `4096 x 256` 的 `uint16` transition table，单这一项约 2 MiB。论文的 FPGA 资源不能由这个数字推断；论文只公开了总体逻辑资源利用率，并未给出精确 RAM 使用量。
- 当前 **SVDCS packet** 同时保留 `mc` 和最终 `payload` 是为了调试验证；实际边端存储/传输应只保留最终压缩载荷和最小元数据。

## 复现边界

- `DWT-Hybrid` 是包内工程实现，不是某篇 DWT 论文的 bit-identical 复现。
- `MMC` 保留论文/作者代码的多模型 + 率约束两阶段框架；v1 residual coder 仍是固定预算稀疏 DCT/Haar-DWT，不是 Antonini/Khan 位平面算术编码逐行复刻。
- `SVDCS` 的 MAX/MIN/SUM/CNT + 阈值 + MC + LZW 主链路按论文实现；论文没有给出 MC 的确切二进制 packing，因此本包采用透明可逆的自定义序列化。

## 报告

- `report/FiveCodec_Experiment_Report.pdf`

报告仍保持两大章：

1. 五算法理论原理 + 压缩后保留参数 + 计算复杂度/资源结构；
2. 压缩比、重构误差、Rate-Distortion、耗时、RAM/空间等实验对比。

Recommended start:

```matlab
preflight_check
run_all_experiments
```

## C port

A toolbox-free C11 translation of the five codecs (ASBC, DWT-Hybrid, CS-OMP, MMC, SVDCS) lives in `../c`. Build and run the encode/decode verification with:

```bash
cd ../c && make test        # IEEE 1159 synthetic
cd ../c && make test-data   # field CSV in ../data
```

See `../c/README.md`.
