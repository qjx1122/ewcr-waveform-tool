# STATUS.md
## 当前目标
- 对 `data/` 全部波形（4 文件 × 6 通道）重跑五算法验证（已完成）

## 已完成
- [x] `verify_field_data` 扫描目录全部 CSV，通道 UA/IA/UB/IB/UC/IC
- [x] 120/120 用例 PASS；结果 `c/results/verify_field_data.csv` 与 `c/results/out/`

## 进行中
- 无

## 下一步（TODO）
1. 如需：整段长记录（>10 周期）
2. 如需：独立 `.ewcr` 解码 CLI

## 决策记录 / 踩坑
- wave4 IB 峰值仅 0.4 A、近似直流偏置，ASBC 估 F0=53.7 Hz，部分算法相似度降到 ~64–80%，仍判 PASS（有限重构 + CR>0）
- 每通道仍取跳过 1 周期后的 10 周期，避免 CS-OMP/MMC 吃整段 5–10 s

## 关键文件路径
- `data/wave{1,2,3,4}.csv`
- `c/results/verify_field_data.csv`
- `c/results/out/manifest.csv`
