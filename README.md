# boltz2-eval — Boltz-2 蛋白-配体亲和力头的实测评测（DAVIS/Vina  docking 子集）

对 Boltz-2（`boltz==2.0.3`）的 affinity prediction head 在同一批 AutoDock Vina
已打分行上做的**同靶点排序 / MSA 消融 / 采样-质量曲线 / 晶体 RMSD / 协议等价性**
评测。21 次 GPU 实测运行 + 6 对原单序列基线，全部原始输出（mmcif +
confidence/affinity JSON + records）与确定性分析脚本随仓库发布，
`results/ext_results.md` 可在无 GPU 条件下由已发布产物原样再生
（本仓库提交前实测：连跑两次哈希不变；与 2026-11-17 原始记录在行尾符归一后逐字节一致，4 张分析表逐字节一致）。

> 结论基调是**负结果 + 口径更正**，这不是失败，是这个仓库的主要价值：
> 见下文「关键结论」与 `docs/MSA口径更正_20261117.md`。

## 目录结构

```
scripts/            全部代码（make_yaml → make_ext_inputs → fetch_msa → run_ext → analyze_ext）
data/dti_subset/    上游 docking 基准的最小切片（panel/pairs/ligands/3 个晶体 PDB）
data/baseline/      原 6 对运行的记录行（preset=affinity-default+msa，实测为单序列，见口径更正）
inputs/             6 个基线 YAML（msa: empty）+ 21 个本批 YAML（msa: msa/*.a3m 或 0/empty）
msa/                MMseqs2 检索并持久化的 3 个 a3m（运行间复用，保证同构建序列同一 MSA）
outputs/            本批 21 次运行：预测 cif、affinity/confidence JSON、processed records JSON
outputs_baseline/   原 6 次运行：同上子集（ss_ 基线，records 里 msa_id=-1 可查）
results/            ext_runs.csv（逐 run 指标）＋ 4 张分析表 ＋ ext_results.md ＋ batch.log
docs/               2026-10-14 原始报告、MSA 口径更正、three-way 判定
```

## 快速复现

结果层（CPU，仅需 pandas/scipy/numpy/biopython）：

```bash
python scripts/analyze_ext.py      # 重新生成 results/ext_results.md 与 ext_*.csv
```

预测层（需要 CUDA GPU，峰值 ~11 GiB，权重缓存 ~6.5 GB 首次自动下载）：

```bash
pip install -r requirements.txt
python scripts/make_yaml.py         # 基线 6 YAML（msa: empty）
python scripts/make_ext_inputs.py   # 本批 21 YAML
python scripts/fetch_msa.py         # 可选：重新检索 a3m（发布件已带）
python scripts/run_ext.py           # 全部 21 次运行（seed 42，可 --only= 子集，断点续跑）
```

从仓库根目录运行。路径常量可用环境变量覆盖：`DTI_BENCH`（默认 `data/dti_subset`）；
`run_ext.py` 用 `BOLTZ_PY`/`BOLTZ_CACHE` 指定解释器与权重缓存。

## 实验设计

| 组 | 问题 | 设计 |
|---|---|---|
| E1 | 同靶点内 Boltz-2 亲和力能否排序真实 pKd | BRAF(8c7x) 全部 6 个 pKd>5.0 配体 + 5 个 pKd=5.0 且 Vina 最强的"对抗"配体，共 11 ligand × 同一持久化 a3m |
| E2 | 真实 MSA vs 单序列的差 | 原 6 对（更正后＝单序列基线）逐一补跑真实 MSA；Δ=真实−单序列；E2rep 重跑 1 对验证可复现 |
| E3 | 结构采样档-质量曲线 | 25步/1回收、100步/2回收、默认200步/3回收；亲和力头固定默认 5×200 |
| E4 | 预测复合体相对晶体的形变 | Cα Kabsch 对齐后整体与结合位点(≤5 Å)RMSD |
| EQ | 服务器 MSA 与本地 a3m 复用是否等价 | 同复合物两路径逐指标对比 |

## 关键结论（详见 `results/ext_results.md`）

- **E1 负结果**：11 配体上 Boltz-2 亲和力 vs pKd Spearman **−0.543**（95% CI [−0.857, +0.080]），Vina −0.286 [−0.812, +0.411]——两者区间都跨 0，同靶点排序在小样本 n 下均未建立。绝对值与 pKd 尺度不可比，只评排序。
- **EQ 方法结论**：本地持久化 a3m 与 `--use_msa_server` 路径结果**逐位相同**（affinity Δ +0.0000，ptm/iptm 相同），且省 1845 s/对 ⇒ 多配体评测应复用 MSA。
- **E2**：加入真实 MSA 后置信度指标有正有负（d_ptm −0.10~+0.09），affinity 绝对值移动可达 1.7；显存约 6100→11000 MiB，单对耗时增至默认档的 2–7 倍；E2rep 重跑逐位复现（seed 固定可复现）。
- **E3**：ptm/iptm 随 steps×recycles 单调上升（如 T0004_L01091 0.419→0.566→0.672）。
- **E4 负结果**：对晶体 Cα RMSD 全部组合整体 18.5–23.7 Å，结合位点 12.7–24.1 Å——短构建序列、无模板输入下 Boltz-2 复合体姿态与晶体差很大；真实 MSA 档的 CA pLDDT 均值约 33–46，不足以自洽预警这种偏差。
- **口径更正（本批起因）**：原 6 对 YAML 写了 `msa: empty`，该字段解析为 `msa_id=-1` 且**优先于** `--use_msa_server`（schema.py L1118 / main.py L517-529），records 证实——原以为是真实 MSA 的运行实为单序列。发布件里 `outputs_baseline/*/*/processed/records/` 可自查。

计时口径：同协议同输入墙钟在夜间长负载下波动 236–4712 s（非协议差异，MSA 张量形状逐位一致），故本仓库只报区间不报单点速度，见 `results/ext_results.md` 末节。

## 数据来源与第三方声明

- 蛋白序列与 pKd/Vina 真值取自 `dti-affinity-bench` 的 docking 结果切片（DAVIS 靶点/配体 + AutoDock Vina 打分 + RCSB PDB 晶体 8c7x/7ukz/6yz4）；未重新对接、未重新取数。
- MSA 由 MMseqs2（boltz 内置协议，经 api.colabfold.com 检索）产生。
- 模型为 MIT 许可的 `boltz==2.0.3` 官方权重；本仓库不含权重。
- 本仓库评测的是**公开模型的可用性**，模型本身非作者开发。

## 未随仓库发布（及原因）

- 模型权重缓存（~6.5 GB，首次运行自动下载）；
- 逐 run 的 pae/pde/plddt npz、processed 缓存 npz/pkl、lightning_logs（可由 `run_ext.py` 重跑再生，指标 JSON 已含结论所需数值）；
- `child.log`（内含本机绝对路径）；
- pKd=5.0 截断组的其余 9 个配体（未跑）。

## 相关仓库

- [esm2-flu-s](https://github.com/fakenice/esm2-flu-s)：同一作者在蛋白语言模型方向的可复现实验仓库。

## License

MIT（见 `LICENSE`）。上游数据的版权归原作者与 RCSB PDB / DAVIS。
