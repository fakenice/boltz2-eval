# Changelog

本仓库演进的如实记录。条目日期 = git 提交日期或所引文档自身的日期；`git log` 可逐条对上。

## [Unreleased]

- 2026-11-20 补齐 `pKd = 5.0` 截断组其余 9 个配体（此前点名的未跑部分），E1 由 11 配体增至 20，
  GPU 实测 9 跑（同协议：本地持久化 a3m、默认结构档、亲和力头默认 5×200、seed 42；
  耗时 226–543 s，峰值显存 10893–11015 MiB，9/9 全部产出）。结果无论正反并入同一张表，并如实记下
  判读变化与三条限制（详见 `docs/DECISIONS.md` 第 9–12 条）：
  **① 方向口径核定**——`affinity_pred_value` 是 log10(IC50)（μM），官方原文 “Lower values indicate
  stronger predicted binding”（jwohlwend/boltz `docs/prediction.md`「Output」节，2026-11-20 取），
  故 ρ(affinity, pKd) 为负才是模型在工作的方向；此前只判「区间是否跨 0」、未对符号下结论。
  **② n=20 全部 E1**：亲和力 ρ=−0.537，95% CI [−0.768, −0.072]（不跨 0），Mantel 置换双侧 p=0.0171，
  随机排序 5% 分位线 |ρ|≥0.447；Vina +0.336 [−0.083, +0.677] p=0.1455、二元结合概率 +0.362
  [−0.088, +0.693] p=0.1191 仍不显著。
  **③ n=11（已发布首批）逐字未变**：亲和力 −0.543 [−0.857, +0.080] p=0.0848，Vina −0.286 p=0.3853，
  随机基线 0.610 ⇒ 与随机排序不可区分。两轮结论不同是**被比较的配体集合**变了（首批截断组只留
  Vina 最强的 5 个），不是模型行为变了；且 20 个 pKd 中 14 个并列在 5.0 上限，n=6 仅活性子集
  什么都判不出（亲和力 −0.257 p=0.6330），3 分数 × 3 子集共 9 个统计量未做多重比较校正。
  代码：`scripts/e1_subsets.py` 单点定义三个子集（名字掩码，不重排行，保证 p 值可复现）；
  `analyze_ext.py`/`permutation_test.py` 子集扩展并对已发布行加漂移断言；
  `check_posebusters.py` 覆盖对象随之由 27 增至 36（B2 结论不变：36/36 最近原子对 0.02–0.27 Å、0/36 建图成功），
  并把写死的「27 个预测」「21 个本批预测」改为按实际行数生成。
  产物：`ext_runs/ext_meta/ext_rmsd/ext_ranking/posebusters_geometry` 五个 csv 只追加、
  原始行按行文本比对逐字保留；`permutation.csv` 增 n=20 一档、子集名 `all_11`→`published_11`；
  `ext_results.md` 的 E1 表加 `batch` 列并新增方向口径段与补齐批计时段。
  文档：README E1/B2/B3/E4/计时/未发布清单同步更新，并把 B2 与 E1 之间改为「未解释的矛盾」表述；
  `docs/DECISIONS.md` 加第 9–12 条、第 4 条计时区间补白天批
- 2026-11-18 新增 `docs/external_baselines.md`：Boltz-2 论文（Pearson 0.66/3-of-8/0.83）与
  PoseBusters 原文数字逐条核对＋口径差异清单；Blaese/J-DTA 未核到、明确不引用；README 加「外部对照」节
- 2026-11-18 新增 `docs/DECISIONS.md`：MSA 口径更正、计时只报区间、a3m 复用、
  B2 有效性/准确度分层、B3 置换协议等决定的仓内依据逐条登记（日期对齐 git）

## 前史（文档留痕，先于本仓库 git 历史）

- 2026-10-14 `docs/报告_本机蛋白结构预测实测_Boltz2_20261014.md`：本机首次实测评测报告
  （21 次 GPU 运行中的早期批次；单序列基线口径的成本与失败率即出自该记录）
- 2026-11-17 `docs/MSA口径更正_20261117.md`：MSA 协议更正——原 6 对 YAML 实为
  `msa: empty`（`msa_id = −1` 优先于 `--use_msa_server`），据此把「6 对全部真实 MSA」
  的早期措辞改为「单序列基线口径」，成本数字只报区间

## [v1.0] - 2026-11-18

仓库建立并当日冻结为可引用版本（Release v1.0，含 tag 与 `.zenodo.json`）。

### Added

- `f12779e` 建仓：21 次 GPU 实测运行 + 6 对单序列基线全部原始输出（mmcif +
  confidence/affinity JSON + records）与确定性 CPU 分析脚本入库
- `776aa71` B2：27 个已发布姿态的 PoseBusters 无参考装载检查（0/27 通过，含正对照原子距）
- `dd1b44a` B3：Mantel 置换检验（20,000 次）+ 随机排序基线
- `97c224f` `CITATION.cff` + README「引用（How to cite）」节
- `cc0883b` `.zenodo.json` 发布元数据

### Changed

- `0accc0d` / `614d155` 报告脱敏：移除本机绝对路径；`.gitattributes` 固定 LF，
  保证跨行尾符逐字节可比
