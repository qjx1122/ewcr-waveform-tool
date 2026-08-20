# STATUS.md
## 当前目标
- 输出「波形压缩比」与「波形相似度」（原波 vs 复原波）（已完成本会话）

## 已完成
- [x] `EwcrMetrics` 增加 Pearson `corr` 与 `similarity_pct`（原波 vs `xhat`）
- [x] 终端/CSV 同时给出压缩比 `CR`/`compress_ratio` 与相似度 `%`
- [x] 合成 25/25、现场 40/40 验证仍全部 PASS

## 进行中
- 无

## 下一步（TODO）
1. 如需：把原波/复原波逐点 CSV 也落盘
2. 如需：全通道或更长现场片段

## 决策记录 / 踩坑
- 压缩比仍用 `CR = (N×16)/bits_compressed`（相对 16 bit PCM）
- 波形相似度用 Pearson `r(x,xhat)`，百分数 `similarity_pct = 100*r`；SNR_dB 保留为能量保真度，不再当作「相似度」
- 仿射变换（增益+偏置）Pearson=1，但 SNR 会下降，两者互补

## 关键文件路径
- `c/src/ewcr_math.c` — `ewcr_unified_metrics`
- `c/results/verify_five_codecs.csv`
- `c/results/verify_field_data.csv`
