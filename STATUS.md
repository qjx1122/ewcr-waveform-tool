# STATUS.md
## 当前目标
- 按 `BOOTSTRAP.md` 完成开局/收尾仪式，固化会话上下文
- 等待用户下达 `ewcr-waveform-tool`（五算法电气波形压缩 MATLAB 基准）的具体任务

## 已完成
- [x] 拉取仓库现状：`arena/01a01e0b-ewcr-waveform-tool` 基于 `b4fc4a3`（首次提供 matlab 原始代码及文档）
- [x] 从 `origin/main`（`3400898`）同步 `BOOTSTRAP.md` 到本工作分支
- [x] 创建 `STATUS.md` 骨架并写入开局状态
- [x] GitHub 认证可用（`gh` + https）
- [ ] MATLAB 实验复现（本环境无 MATLAB / Octave）

## 进行中
- 会话引导协议落地（开局 + 收尾）

## 下一步（TODO）
1. 等用户明确本会话任务（算法改动 / 复现实验 / 文档 / 移植等）
2. 若需跑基准：准备 MATLAB R2025b + Signal/Wavelet/Communications Toolbox，或评估 Octave 替代
3. 对照 `matlab/README.md` 的一键流程：`preflight_check` → `run_all_experiments`

## 决策记录 / 踩坑
- 本工作分支由 Arena session 固定为 `arena/01a01e0b-ewcr-waveform-tool`，不按 BOOTSTRAP 另开 `feature/xxx`，以免会话丢失
- `origin/main` 比本分支多 1 个提交：`3400898 上传BOOTSTRAP.md文件`；已把该文件检出到本分支，未整分支 rebase/merge，避免改写 Arena 跟踪分支
- 仓库无 `setup.sh` / Node 依赖；工程主体在 `matlab/`，推荐环境是 MATLAB R2025b + 三个 Toolbox
- 本沙箱未安装 `matlab` / `octave`，无法执行 `preflight_check` / `run_all_experiments`
- 会话纪要只追加到 `session/NILM_AC_session_complete.md`，专题报告只追加到 `REPORT_TEST.md`，不另建文件

## 关键文件路径
- `BOOTSTRAP.md` — 开局/收尾协议
- `STATUS.md` — 本文件
- `session/NILM_AC_session_complete.md` — 会话纪要（只追加）
- `matlab/README.md` — MATLAB 基准说明与运行方式
- `matlab/run_all_experiments.m` — 一键入口
- `matlab/preflight_check.m` — Toolbox 预检
- `matlab/algorithms/{ASBC,DWT_Hybrid,CS_OMP,MMC,SVDCS}/` — 五算法
- `matlab/experiments/` — exp1–exp6
- `matlab/results/` — CSV/MAT 结果（含 reference）
- `matlab/figures/` — 实验图
- `matlab/report/FiveCodec_Experiment_Report.pdf` — 实验报告
