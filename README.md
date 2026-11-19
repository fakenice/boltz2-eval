# boltz2-eval — Boltz-2 蛋白-配体亲和力头的实测评测（DAVIS/Vina  docking 子集）

对 Boltz-2（`boltz==2.0.3`）的 affinity prediction head 在同一批 AutoDock Vina
已打分行上做的**同靶点排序 / MSA 消融 / 采样-质量曲线 / 晶体 RMSD / 协议等价性**
评测。30 次 GPU 实测运行（21 次为 2026-11-17/18 首批，9 次为 2026-11-20 补齐
pKd=5.0 截断组）+ 6 对原单序列基线，全部原始输出（mmcif +
confidence/affinity JSON + records）与确定性分析脚本随仓库发布，
`results/ext_results.md` 可在无 GPU 条件下由已发布产物原样再生
（本仓库提交前实测：`analyze_ext.py`、`permutation_test.py`、`check_posebusters.py` 各自连跑两次产物哈希不变。
2026-11-20 补齐对已发布产物是**只追加**：`ext_runs/ext_meta/ext_rmsd/ext_ranking/posebusters_geometry` 五个 csv 的 2026-11-17/18 原始行按行文本比对全部逐字保留、表头未变，
`permutation.csv` 里 n=11 与 n=6 两档数值不变（只把子集名 `all_11` 改为 `published_11` 并新增 n=20 一档）。）

> 结论基调是**负结果 + 口径更正**，这不是失败，是这个仓库的主要价值：
> 见下文「关键结论」与 `docs/MSA口径更正_20261117.md`。
> 2026-11-20 的补齐让 E1/B3 的判读发生了一处方向变化（同一批分数、被比较的配体集合变了），
> 已在下面如实标注，v1.0 发布件里的 n=11 数字保持原样、可继续引用。

## 目录结构

```
scripts/            全部代码（make_yaml → make_ext_inputs → fetch_msa → run_ext → analyze_ext；check_posebusters 为姿态有效性检查；permutation_test 为置换检验＋随机基线）
data/dti_subset/    上游 docking 基准的最小切片（panel/pairs/ligands/3 个晶体 PDB）
data/baseline/      原 6 对运行的记录行（preset=affinity-default+msa，实测为单序列，见口径更正）
inputs/             6 个基线 YAML（msa: empty）+ 30 个本批 YAML（msa: msa/*.a3m 或 0/empty）
msa/                MMseqs2 检索并持久化的 3 个 a3m（运行间复用，保证同构建序列同一 MSA）
outputs/            本批 30 次运行：预测 cif、affinity/confidence JSON、processed records JSON
outputs_baseline/   原 6 次运行：同上子集（ss_ 基线，records 里 msa_id=-1 可查）
results/            ext_runs.csv（逐 run 指标）＋ 4 张分析表 ＋ ext_results.md ＋ posebusters_geometry.csv/md ＋ permutation.csv/md ＋ batch.log（首批与 2026-11-20 补齐批两段，只追加）
docs/               2026-10-14 原始报告、MSA 口径更正、three-way 判定
```

## 快速复现

结果层（CPU，仅需 pandas/scipy/numpy/biopython）：

```bash
python scripts/analyze_ext.py      # 重新生成 results/ext_results.md 与 ext_*.csv
python scripts/permutation_test.py # 重新生成 results/permutation.csv 与 permutation.md（E1 的置换检验＋随机基线）
```

有效性层（CPU，独立小环境，`pip install -r requirements-posebusters.txt`）：

```bash
python scripts/check_posebusters.py  # 重新生成 results/posebusters_geometry.csv 与 posebusters.md
```

预测层（需要 CUDA GPU，峰值 ~11 GiB，权重缓存 ~6.5 GB 首次自动下载）：

```bash
pip install -r requirements.txt
python scripts/make_yaml.py         # 基线 6 YAML（msa: empty）
python scripts/make_ext_inputs.py   # 本批 30 YAML
python scripts/fetch_msa.py         # 可选：重新检索 a3m（发布件已带）
python scripts/run_ext.py           # 全部 30 次运行（seed 42，可 --only= 子集，断点续跑）
```

从仓库根目录运行。路径常量可用环境变量覆盖：`DTI_BENCH`（默认 `data/dti_subset`）；
`run_ext.py` 用 `BOLTZ_PY`/`BOLTZ_CACHE` 指定解释器与权重缓存。

## 实验设计

| 组 | 问题 | 设计 |
|---|---|---|
| E1 | 同靶点内 Boltz-2 亲和力能否排序真实 pKd | BRAF(8c7x) 全部 6 个 pKd>5.0 配体 + pKd=5.0 截断组（首批只取 Vina 打分最强的 5 个，即对亲和力头最不利的挑法；2026-11-20 补齐其余 9 个 ⇒ 截断组整组 14 个），共 20 ligand × 同一持久化 a3m |
| E2 | 真实 MSA vs 单序列的差 | 原 6 对（更正后＝单序列基线）逐一补跑真实 MSA；Δ=真实−单序列；E2rep 重跑 1 对验证可复现 |
| E3 | 结构采样档-质量曲线 | 25步/1回收、100步/2回收、默认200步/3回收；亲和力头固定默认 5×200 |
| E4 | 预测复合体相对晶体的形变 | Cα Kabsch 对齐后整体与结合位点(≤5 Å)RMSD |
| EQ | 服务器 MSA 与本地 a3m 复用是否等价 | 同复合物两路径逐指标对比 |
| B2 | 预测配体姿态是否物理有效 | PoseBusters 口径两层：原始几何指标（最近原子对/vdW 重叠，无键级假设）＋标准装载/建图检查；共晶配体作正对照 |
| B3 | E1 排序与随机基线能否区分 | Mantel 置换（整排 pKd、保留截断并列）双侧 p ＋零分布 |ρ| 95 分位作随机基线；三档子集：全部 E1（n=20）、首批已发布（n=11）、仅活性（n=6），子集定义集中在 `scripts/e1_subsets.py` |

## 关键结论（详见 `results/ext_results.md`）

- **E1（截断组补齐后判读变了一处，详见 `results/ext_results.md`）**：补齐 pKd=5.0 组其余 9 个配体后，20 个 BRAF 配体上 Boltz-2 亲和力 vs pKd Spearman **−0.537**（95% CI [−0.768, −0.072]）——方向口径 2026-11-20 依官方文档核定：`affinity_pred_value` 是 log10(IC50)（μM），官方原文 “Lower values indicate stronger predicted binding”，而 pKd 越大结合越强，**所以负相关才是模型在工作的方向**；这个区间不再跨 0。Vina **+0.336** [−0.083, +0.677]、二元结合概率 +0.362 [−0.088, +0.693]（二者正＝正确方向）区间仍跨 0。同时必须如实记下三条限制：
  - 首批已发布的 n=11 视图（亲和力 −0.543 [−0.857, +0.080]、Vina −0.286 [−0.812, +0.411]）字节不变、可继续引用。两轮结论不同**不是模型行为变了，而是被比较的配体集合变了**：n=11 的截断组只留了 Vina 最强的 5 个（对亲和力最不利的挑法），n=20 是整组 14 个。
  - 20 个 pKd 里有 14 个并列在 DAVIS 的 5.0 截断上限，这个统计量主要衡量"活性配体与截断配体能否分开"，不是活性配体之间的精细排序；仅活性子集（n=6）亲和力 −0.257 [−1.000, +1.000]，什么也判不出。
  - 共考察了 9 个统计量（3 个分数 × 3 个子集）且未做多重比较校正，p=0.0171 经不起这个量级的挑选；本仓库把它记为**方向上与真实 pKd 一致、但远不足以作为可用性结论**。绝对值与 pKd 尺度不可比，只评排序。
- **EQ 方法结论**：本地持久化 a3m 与 `--use_msa_server` 路径结果**逐位相同**（affinity Δ +0.0000，ptm/iptm 相同），且省 1845 s/对 ⇒ 多配体评测应复用 MSA。
- **E2**：加入真实 MSA 后置信度指标有正有负（d_ptm −0.10~+0.09），affinity 绝对值移动可达 1.7；显存约 6100→11000 MiB，单对耗时增至默认档的 2–7 倍；E2rep 重跑逐位复现（seed 固定可复现）。
- **E3**：ptm/iptm 随 steps×recycles 单调上升（如 T0004_L01091 0.419→0.566→0.672）。
- **E4 负结果**：30 个预测（24 真实 MSA + 6 单序列）对晶体 Cα RMSD 全部组合整体 18.5–23.8 Å，结合位点 12.7–24.1 Å——短构建序列、无模板输入下 Boltz-2 复合体姿态与晶体差很大；真实 MSA 档的 CA pLDDT 均值约 33–46，不足以自洽预警这种偏差。
- **B2 负结果（姿态有效性）**：36/36 个预测（2026-11-20 起含补齐的 9 个；此前 27/27 结论相同）的配体原子云中都存在 0.02–0.27 Å 的原子对——低于任何真实共价键一个量级；共晶配体正对照同一指标 1.14–1.27 Å，排除解析伪影。按 PoseBusters 协议从坐标感知连通图并做模板键级赋值，0/36 能建成合法分子图（全部卡在装载层）。明细 `results/posebusters.md`。这与补齐后 E1 的方向一致**不互相印证，而是一处未被解释的矛盾**：同一批运行里，姿态层面不物理（配体原子重叠、建不出合法分子图、对晶体 Cα RMSD 18–24 Å），打分层面在 n=20 上却与真实 pKd 的排序方向一致（其限制见 E1 条）。两条结论各自都由已发布产物支撑，本仓库只如实并列，不给机制解释。
- **B3（统计口径补强，与 E1 同步更新）**：Mantel 置换（整排 pKd、保留截断并列，20000 次、seed 42、+1 校正）三档子集结果不同：n=20 时 Boltz-2 亲和力 ρ=−0.537 双侧 p=0.0171，随机排序有 5% 概率碰出 |ρ|≥0.447，观测值在这条线之上；Vina ρ=+0.336 p=0.1455、二元结合概率 ρ=+0.362 p=0.1191 均不过 0.05。n=11（已发布首批）时三者全不显著——亲和力 ρ=−0.543 p=0.0848、Vina ρ=−0.286 p=0.3853、二元 ρ=+0.353 p=0.2715，且随机 5% 分位线是 |ρ|≥0.610，三个观测值都在线下 ⇒ 与随机排序不可区分。n=6 仅活性子集只能给粗判（亲和力 p=0.6330）。
  n=11 与 n=20 的差别同样是配体集合变化所致，不是同一假设上的两次独立检验；叠加上面 9 个统计量未做多重比较校正，不要把 n=20 的 p=0.0171 单独当作"Boltz-2 排序显著"的证据引用。两个方向都不构成「显著差于随机」的证据。详见 `results/permutation.md`。
- **口径更正（本批起因）**：原 6 对 YAML 写了 `msa: empty`，该字段解析为 `msa_id=-1` 且**优先于** `--use_msa_server`（schema.py L1118 / main.py L517-529），records 证实——原以为是真实 MSA 的运行实为单序列。发布件里 `outputs_baseline/*/*/processed/records/` 可自查。

计时口径：同协议同输入的单对墙钟在全部运行间波动 226–4712 s，与计算内容无关（MSA 张量形状逐位一致，EQ 与 E1 同复合物结果逐位相同）；其中 2026-11-17/18 夜间批为 236–4712 s，2026-11-20 补齐批在白天独占显卡下是 226–543 s，全部落在夜间批上界之下，与「慢档来自长时间负载下的 GPU 降频或后台进程争用」一致（未单独确证，是相关而非因果）。故本仓库只报区间不报单点速度，见 `results/ext_results.md` 末节。

## 外部对照（已发表数字，核过原文）

完整对照与逐条口径差异见 [`docs/external_baselines.md`](docs/external_baselines.md)
（2026-11-18 逐条抓取原文核对）。要点：

- **Boltz-2 论文全文不报 Spearman**（检索 0 次命中）：亲和力指标是 per-assay Pearson——
  FEP+ 4-target 子集均值 0.66、Recursion 8 盲测中 3 个 > 0.55、TYK2（Hahn 基准）R = 0.83。
  与本仓库 within-target Spearman **统计量、集合、粒度都不同**，只作方向参照，不并表。
- **PoseBusters 原文同样零亲和力指标**，只报 RMSD 与有效性比例——支撑本仓库 B2 把
  "有效性"与"准确度"分开的协议设计；不存在可引的"Vina 在 PoseBusters 上的 Spearman"。
- **DAVIS within-target Spearman 无已发表基线可引**：常被提及的 Blaese et al. 2022 与 J-DTA
  两条一手全文均未核到（检索路线与失败事实记录在文档 §4），本仓库不引用任何未核数字。

## 数据来源与第三方声明

- 蛋白序列与 pKd/Vina 真值取自 `dti-affinity-bench` 的 docking 结果切片（DAVIS 靶点/配体 + AutoDock Vina 打分 + RCSB PDB 晶体 8c7x/7ukz/6yz4）；未重新对接、未重新取数。
- MSA 由 MMseqs2（boltz 内置协议，经 api.colabfold.com 检索）产生。
- 模型为 MIT 许可的 `boltz==2.0.3` 官方权重；本仓库不含权重。
- 本仓库评测的是**公开模型的可用性**，模型本身非作者开发。

## 未随仓库发布（及原因）

- 模型权重缓存（~6.5 GB，首次运行自动下载）；
- 逐 run 的 pae/pde/plddt npz、processed 缓存 npz/pkl、lightning_logs（可由 `run_ext.py` 重跑再生，指标 JSON 已含结论所需数值）；
- `child.log`（内含本机绝对路径）。

v1.0 发布件里这里还列着「pKd=5.0 截断组的其余 9 个配体（未跑）」——2026-11-20 已跑完并入 E1，本仓库不再有未跑项。

## 相关仓库

- [esm2-flu-s](https://github.com/fakenice/esm2-flu-s)：同一作者在蛋白语言模型方向的可复现实验仓库。

## 引用（How to cite）

本仓库以 GitHub Release v1.0 作为固定引用版本；机器可读元数据见根目录 `CITATION.cff`（GitHub 会自动在侧栏展示引用建议）。引用负结果与引用正结果同样重要。若后续分配了 DOI，请优先引用 DOI。

```bibtex
@software{wu2026boltz2eval,
  author  = {{吴奕婷 (Yiting Wu)}},
  title   = {boltz2-eval：Boltz-2 蛋白-配体亲和力头的单卡严格评测（DAVIS / AutoDock Vina 子集）},
  url     = {https://github.com/fakenice/boltz2-eval},
  version = {1.0},
  date    = {2026-11-18}
}
```

## License

MIT（见 `LICENSE`）。上游数据的版权归原作者与 RCSB PDB / DAVIS。
