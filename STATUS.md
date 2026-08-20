# STATUS.md
## 当前目标
- 输出压缩文件与复原波形文件（已完成）

## 已完成
- [x] 五算法打包 `.ewcr` 压缩码流
- [x] 写出 `*_original.csv` / `*_reconstructed.csv` / `*_compare.csv`
- [x] 现场验证 `results/out/manifest.csv` 索引
- [x] 合成 25/25、现场 40/40 仍全部 PASS

## 进行中
- 无

## 下一步（TODO）
1. 如需：从 `.ewcr` 独立解码的命令行工具
2. 如需：全通道 / 更长现场片段

## 决策记录 / 踩坑
- 压缩文件为自定义小端容器，魔数 `EWCR`；CSV 为可读波形
- `results/out/` 下大批 CSV/bin 不进 Git（gitignore），只提交 README 与 manifest
- 验证仍在内存编解码；落盘是额外产物，不改变算法

## 关键文件路径
- `c/results/out/` — 压缩与复原文件
- `c/src/ewcr_io.c` — 读写
- `c/results/out/manifest.csv`
