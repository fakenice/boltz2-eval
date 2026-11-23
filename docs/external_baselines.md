# 外部对照（已发表数字，逐条核过原文）

核验日期 2026-11-18；取证方式：抓取论文原文全文（Europe PMC fullTextXML / NCBI PMC）后
对关键句与表号做字符串级核对。本仓库与发表工作**没有一项同口径的可比数字**——
这本身就是本文件要如实说明的事：已发表文献里没有"DAVIS 单靶点、小配体集、within-target
Spearman"这一格的已发表数字（不可核到的条目见 §4，一律不引用）。

## 1. 本仓库锚（先固定自己的口径）

- 指标：**within-target Spearman**（BRAF 单靶点，pKd 含 5.0 截断并列），bootstrap CI。
- 代表数字（2026-11-18 首版口径，n=11 ＝ 首批 6 活性 + Vina 最强 5 截断）：
  Boltz-2 亲和力 ρ = **−0.543** [−0.857, +0.080]；Vina ρ = −0.286；
  Mantel 置换双侧 p = 0.0864（原 0.0848，2026-11-20 修复等值置换漏计）；姿态有效性 0/27 建图通过（B2）；
  对晶体 Cα RMSD 整体 18.5–23.7 Å（E4）。
- 2026-11-20 补齐截断组其余 9 个配体后的同一口径（n=20）：亲和力 ρ = **−0.537** [−0.768, −0.072]，
  置换 p = 0.0171；Vina +0.336 p = 0.1465（原 0.1455，同日计数修复）；姿态有效性 0/36；E4 整体 18.5–23.8 Å。
  方向口径同日核定：`affinity_pred_value` 是 log10(IC50)（μM），越小结合越强
  （jwohlwend/boltz `docs/prediction.md`「Output」节）⇒ **ρ 为负才是模型在工作的方向**。
  引用本文件时请注意：n=11 与 n=20 的区别是被比较的配体集合不同（n=11 的截断组是被挑过的一半），
  不是模型行为变了；两者的限制与多重比较问题见 `docs/DECISIONS.md` 第 9–10 条。
- 配置：`boltz==2.0.3`、seed 42、affinity 头默认 5×200、E1 用单一持久化 a3m。

## 2. Boltz-2 论文（Passaro, Corso, Wohlwend et al., bioRxiv 2025.06.14.659707）

出处锚：bioRxiv 预印本（无 arXiv 版）；经 Europe PMC PPR1039145 全文核对，下列句段均逐字在文中。

- **全文 "Spearman" 0 次命中**——论文不报告任何 Spearman，亲和力指标是
  "Pearson correlation averaged over each assay"（Fig. 6 及其 legend）。
- 正文给出的亲和力数字（逐条核对）：
  - FEP+ 4-target 子集：**平均 Pearson 0.66**（"On the 4-target FEP subset, Boltz-2
    achieves an average Pearson correlation of 0.66"）。
  - Recursion 8 个盲测 assay：**3/8 达到 Pearson > 0.55**，其余 5 个"有限"
    （原文并指出这种 assay 间波动在 FEP 方法里也常见，GPCR 类靶点尤其弱）。
  - TYK2：基于 Hahn et al. 2022 基准数据 **Pearson R = 0.83**；
    Boltz-ABFE（显式 FEP 管线，非 Boltz-2 本体）R = 0.95、MAE 0.42 kcal/mol。
- 评测集：四个 affinity test sets = OpenFE 876 hit-to-lead、FEP+ 4-target 子集、
  CASP16 affinity（2 靶点 140 对）、Recursion 8 assay。**不含 DAVIS within-target 排序基准。**
  MSA/steps/repeat 等推理超参正文未给出（未核到，不外推）。
- 口径差异（引用上表数字前必读）：
  1. 统计量不同：per-assay Pearson（对 log10 µM 绝对值线性相关；与本仓库 2026-11-20 核定的
     `affinity_pred_value` = log10(IC50/μM) 同一量纲，故论文与本仓库评的是同一族量的不同统计量）
     vs 本仓库 within-target
     Spearman（只看排序，尺度无关）。
  2. 集合不同：4 组 assay 均值/加权 vs 本仓库单靶点 11 配体（2026-11-20 起 20 配体，pKd 截断并列）。
  3. 协议不同：论文 out-of-the-box 全流程 vs 本仓库固定 seed、复用 a3m 的评测口径。
- 方向性判读：论文自报"8 个盲测 assay 只有 3 个 Pearson > 0.55"说明 affinity 头的
  assay 间方差本来就大；本仓库在单靶点截断小样本上，n=11 未检出关联（不等于证明随机或等效），补齐截断组后的
  n=20 亲和力排序与其同向且区间不跨 0（限制见 §1 与 `docs/DECISIONS.md` 第 10 条），
  与该方差叙述**都不矛盾**——方差大既能容下"某一半样本碰不出信号"，也能容下"整组样本出现同向排序"，
  但既不构成对论文结果的复现，也不构成反驳。

## 3. PoseBusters 论文（Buterez et al., Chemical Science 2024；arXiv:2308.05777 / PMC10901501）

经 PMC 全文核对：**Spearman / Pearson / Kendall / affinity 均 0 次命中**——
该论文只报姿态 RMSD ≤ 2 Å 比例与 PB-Valid 通过比例（AstexDiverse 85 例 +
PoseBusters Benchmark 集），不做亲和力评测。因此不存在"Vina 在 PoseBusters 上的
within-target Spearman"可引。

与本仓库 B2 的真实关系是**方法论支撑**：论文核心论点是"RMSD 小不等于姿态物理有效"，
有效性检查必须独立做；本仓库 B2 因此只判"能否从坐标建成合法分子图"，
0/27（2026-11-20 起 0/36）卡在装载层正是该论文警告的失效类型在 Boltz-2 预测姿态上的体现。

## 4. 未核到项（本仓库一律不引用其数字）

- **Blaese et al. 2022《Affinity comparisons: are we there yet?》（MLComp）**：该方向常被引作
  DAVIS within-target Spearman 的基线参照，但 2026-11-18 未能获取一手全文
  （arXiv API、Europe PMC 标题检索、Crossref、chemRxiv 各路线均失败或未命中此文；
  网传 arXiv:2111.02638 不是该文）。**没有核到原文前，本仓库不写它的任何数字。**
- **J-DTA（被声称的 ICML 2024 联合 docking+affinity 工作）**：arXiv 全文检索与
  ICML 2024 接收列表均未命中，出处未核到，不引用。
- 易错锚点声明：PoseBusters 正身是 arXiv:2308.05777（2310.02312 非此文）；
  Boltz-2 无 arXiv 版，引用用 bioRxiv DOI；合作者列表以 bioRxiv 页面为准。

## 5. 这组对照支撑什么、不支撑什么

- 支撑：把"发表数字与本仓库数字不可并列"写成显式口径差异，而不是硬凑一张对比表。
- 支撑：B2 的协议选择（ validity 与 accuracy 分离）有 PoseBusters 原文依据。
- 不支撑：用 Boltz-2 论文的 Pearson 0.66/0.83 说"本仓库 −0.543 与文献矛盾"或
  "文献显示该头很强"——统计量、集合、粒度都不同。
- 不支撑：由本仓库结果外推"Boltz-2 affinity 头在 DAVIS 上失效"——只测了 1 靶点 ×
  11 配体（2026-11-20 起 20 配体）、一种固定配置；负结果的边界写在 README「关键结论」各条与 `docs/DECISIONS.md`。
- 同样不支撑（2026-11-20 新增）：把 n=20 的"方向与真实 pKd 一致、p=0.0171"读成"该头可用"——
  14/20 的 pKd 并列在 5.0 截断上限、9 个统计量未做多重比较校正、且同一批姿态在 B2/E4 上仍不物理
  （`docs/DECISIONS.md` 第 9–10 条；README §关键结论 E1/B2 条）。

## 6. 引用锚（均于 2026-11-18 实际抓取核对）

- Boltz-2：https://www.biorxiv.org/content/10.1101/2025.06.14.659707 （全文经
  Europe PMC PPR1039145 fullTextXML 核对，Fig. 6 caption 与 Results "Affinity" 段逐字核）。
- PoseBusters：https://arxiv.org/abs/2308.05777 ；https://pmc.ncbi.nlm.nih.gov/articles/PMC10901501/
  （fullTextXML 指标词命中计数核）。
- 未核到项：见 §4，只记录检索路线与失败事实，不给数字。
