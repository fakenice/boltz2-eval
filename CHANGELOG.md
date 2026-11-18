# Changelog

本仓库演进的如实记录。条目日期 = git 提交日期或所引文档自身的日期；`git log` 可逐条对上。

## [Unreleased]

- 计划：补齐 `pKd = 5.0` 截断组其余 9 个配体（README「未随仓库发布」一节点名的未跑部分），
  结果无论正反并入同一张负结果表
- 计划：外部对照表（Boltz-2 论文 / Affinity Benchmark / PoseBusters 原文数字，逐条带出处与口径差异）

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
