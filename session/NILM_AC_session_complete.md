# NILM_AC 会话纪要

## [2026-08-20] 会话纪要
- 目标：读取 `BOOTSTRAP.md` 并执行开局仪式与收尾仪式，建立可续接的会话上下文
- 完成项：
  - `git status` / `git log -10 --oneline` / `git branch --show-current`
  - 确认 `STATUS.md`、`session/NILM_AC_session_complete.md`、`REPORT_TEST.md` 原先均不存在；按模板创建前两者
  - GitHub 登录正常；无 `setup.sh`/`pnpm`，跳过 Node 依赖恢复
  - 从 `origin/main` 检出 `BOOTSTRAP.md` 到工作分支
  - 确认本环境无 MATLAB/Octave，无法跑基准
- 关键决策：
  - 保持 Arena 固定分支 `arena/01a01e0b-ewcr-waveform-tool`，不切换/不另建 feature 分支
  - 只检出 `BOOTSTRAP.md`，不 merge 整个 `origin/main`（二者内容等价，仅该文件之差）
  - 尚无用户专题，不创建/不更新 `REPORT_TEST.md` 与 `REPORT.md`
- 未决问题：
  - 用户尚未给出本会话的实质性开发任务
  - 沙箱缺少 MATLAB R2025b 与 Wavelet/Signal/Communications Toolbox
- 相关文件/分支：`BOOTSTRAP.md`、`STATUS.md`、`session/NILM_AC_session_complete.md`、`arena/01a01e0b-ewcr-waveform-tool`

## [2026-08-20] 会话纪要
- 目标：把 `matlab/` 下五算法转换成 C，并用五算法做验证测试
- 完成项：
  - 新增 `c/`：ASBC / DWT-Hybrid / CS-OMP / MMC / SVDCS 编解码 + IEEE 1159 信号 + 验证入口
  - `make test`：8 内核单测 + 4 信号×5 算法（2 周期）+ 纯正弦 4 周期，共 25 例全部 PASS
  - 纯正弦 MMC 与 MATLAB 参考 CR/SNR 对齐；ASBC/CS/SVDCS 高 SNR 重构
  - 写入 `c/results/verify_five_codecs.csv`、`c/README.md`、`REPORT_TEST.md`
- 关键决策：
  - 不依赖 MATLAB/FFTW：自实现 FFT、DCT、Haar、db4、Huffman、LZW
  - 验证是 C 编码器/C 解码器闭环，不要求与 MATLAB bit-identical
  - db4 采用周期延拓正交滤波器组（与 MATLAB `wavedec` 默认 `sym` 不同）
  - CS-OMP 测量矩阵用本地 RNG，CR 公式与 MATLAB 相同（只依赖 M、ybits）
- 未决问题：
  - 未移植 exp1–exp6 全套实验
  - 短窗 DWT/SVDCS 的 CR 受头开销影响可小于 1
- 相关文件/分支：`c/`、`matlab/README.md`、`REPORT_TEST.md`、`arena/01a01e0b-ewcr-waveform-tool`

## [2026-08-20] 会话纪要
- 目标：拉取最新版本，用 `data/` 现场波形对五算法做验证测试
- 完成项：
  - `git fetch` + `reset --hard` 到 `ca8f588 上传现场真实波形数据`
  - 解析 10 kHz 三相 CSV；新增 `ewcr_load_wave_csv` 与 `verify_field_data`
  - `make test-data`：wave1–4 × UA/IA × 五算法，40/40 PASS
  - MMC 在窗长 200 时跳过非 2 幂 Haar 残差
- 关键决策：
  - 每通道取跳过 1 周期后的 10 周期（N=2000），避免 CS-OMP 对整段 5–10 s 爆内存
  - SVDCS 按片段峰值设 `signal_scale`
  - 现场用例不设合成纯正弦的 SNR≥12 dB 门槛，只要求有限重构与 CR>0
- 未决问题：
  - 未跑 UB/IB/UC/IC 与整段长记录
  - ASBC 现场 F0 估计受 FFT 栅格影响（48.83 Hz vs 50 Hz）
- 相关文件/分支：`data/`、`c/src/verify_field_data.c`、`c/results/verify_field_data.csv`、`arena/01a01e0b-ewcr-waveform-tool`

## [2026-08-20] 会话纪要
- 目标：增加输出波形压缩比与波形相似度（原始波形 vs 复原波形）
- 完成项：
  - `ewcr_unified_metrics` 计算 Pearson `corr`、`similarity_pct`
  - 合成/现场验证终端与 CSV 增加 `compress_ratio`、`similarity_pct`、`corr`
  - 单测 corr identity/affine；`make test` 与 `make test-data` 全过
- 关键决策：
  - 压缩比 = CR（16 bit 计账）；相似度 = Pearson，不用 SNR 冒充相似度
  - SNR_dB 仍输出，作为能量保真度
- 未决问题：未落盘逐点 original/reconstructed 波形 CSV
- 相关文件/分支：`c/include/ewcr.h`、`c/src/ewcr_math.c`、`c/src/verify_*.c`、`c/results/*.csv`
