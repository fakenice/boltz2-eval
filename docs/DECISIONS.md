# Decision log

本仓库关键口径与协议决定的记录。每条注明**仓内依据**（文件/代码位置）与可核对的日期；
只收当时已留下书面或代码痕迹的决定，不补写无痕迹的叙事。

## 2026-10-14（建仓前实测期，`docs/报告_本机蛋白结构预测实测_Boltz2_20261014.md` 留痕）

1. **绝对值不可比，只评排序**：Boltz-2 affinity 头输出与 pKd 尺度不同，
   E1 设计定为 within-target Spearman，不并表绝对亲和力
   （原始报告与 `README.md` §关键结论 E1 条）。
2. **单序列基线口径的成本与失败率数字**出自该原始报告的早期批次记录，
   后续引用沿用其口径（`CHANGELOG.md` 前史条）。

## 2026-11-17（`docs/MSA口径更正_20261117.md`）

3. **口径更正：原 6 对"真实 MSA"运行改称"单序列基线"**。证据：YAML 写 `msa: empty`
   解析为 `msa_id=-1`，且该字段**优先于** `--use_msa_server`（boltz `schema.py` L1118 /
   `main.py` L517–529）；records 可自查（发布件 `outputs_baseline/*/*/processed/records/`）。
   这是本批实验的起因，也是"早期措辞不撤回、以更正文档为准"的决定。
4. **计时只报区间、不报单点速度**：同协议同输入墙钟在夜间长负载下波动 236–4712 s
   （非协议差异，MSA 张量形状逐位一致）（`README.md` 计时口径段；`results/ext_results.md` 末节）。

## 2026-11-18（提交 `f12779e`–`97c224f`，建仓日）

5. **多配体评测复用本地持久化 a3m，不用服务器 MSA**。EQ 组实测两路径结果逐位相同
   （affinity Δ +0.0000，ptm/iptm 相同）且省 1845 s/对；同一构建序列同一 MSA 由
   `msa/` 目录保证（`README.md` §关键结论 EQ 条；`inputs/` YAML 指向 `msa/*.a3m`）。
6. **B2 把"有效性"与"准确度"分层检查**：先原始几何指标（最近原子对/vdW 重叠，
   不做键级假设），再 PoseBusters 装载/建图；共晶配体作正对照以排除解析伪影
   （`scripts/check_posebusters.py`；`results/posebusters.md`；
   `docs/external_baselines.md` §3——PoseBusters 原文同样零亲和力指标，支撑此分离）。
7. **B3 置换检验协议**：Mantel 整排置换 pKd、保留截断并列、双侧 p、20,000 次、
   seed 42、+1 校正；并同时报随机基线 |ρ| 95 分位（n=11 时 0.610），
   bootstrap CI 与 p 值双表述（`scripts/permutation_test.py`；`results/permutation.md`）。
8. **发布件脱敏与行尾符**：报告移除本机绝对路径（`0accc0d`/`614d155`）；
   `.gitattributes` 固定 LF，保证复现产物逐字节可比（`CHANGELOG.md` Changed 节；
   `README.md` 首段"逐字节一致"声明的前提）。
