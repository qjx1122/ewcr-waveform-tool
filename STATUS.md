# STATUS.md
## 当前目标
- 拉取最新现场波形，用 `data/` 对五算法做 C 编解码验证（已完成本会话范围）

## 已完成
- [x] 拉取 `ca8f588 上传现场真实波形数据`（`data/wave1..4.csv`）
- [x] C11 五算法 + IEEE 1159 合成验证（此前）
- [x] 现场 CSV 加载器 + `verify_field_data`：4 文件 × UA/IA × 5 算法 = 40 例全部 PASS
- [x] 结果表 `c/results/verify_field_data.csv`

## 进行中
- 无

## 下一步（TODO）
1. 如需：全通道（UB/IB/UC/IC）或更长片段（>10 周期）
2. 如需：把 exp1–exp6 扫频也做到 C
3. ASBC 现场 F0 落在 48.83 Hz（N=2000 FFT 栅格），可改为已知 50 Hz 先验

## 决策记录 / 踩坑
- 现场数据 10 kHz、约 230 Vrms / 1 A，三相 `timestamp,UA,IA,UB,IB,UC,IC`
- 验证片段：跳过首周期后取 10 周期（N=2000）；CS-OMP 不能吃整段 10 s（Φ 过大）
- SVDCS `signal_scale` 取片段峰值，避免 330 V 被 2.0 p.u. 满量程截幅
- MMC 窗长 N=fs/f0=200 非 2 的幂，Haar 残差跳过，只竞 none/DCT
- 现场验证口径：有限重构、CR>0、SNR 有限（不要求合成纯正弦的 12 dB 门槛）

## 关键文件路径
- `data/wave{1,2,3,4}.csv` — 现场波形
- `c/src/verify_field_data.c` — 现场验证入口（`make test-data`）
- `c/results/verify_field_data.csv`
- `c/src/verify_five_codecs.c` — 合成验证（`make test`）
- `c/README.md`
