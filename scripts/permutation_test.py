"""B3: permutation tests and random-ranking baselines for the E1 ranking result.

The E1 subsets in ext_results.md have different bootstrap-CI conclusions.
This script adds a complementary permutation view. Under the null "the score
carries no information about pKd", Spearman rho is exactly distributed over
the pairings of the two rank vectors (Mantel permutation test). Ligands share
the DAVIS truncation value pKd = 5.0, so the pKd vector is permuted as a whole
against the fixed score ranks, which preserves the tie structure.

Monte-Carlo two-sided p-values (20000 shuffles, rng seed 42, +1 correction) for
the three E1 scores (Boltz-2 affinity, -Vina, binary bind probability) on every
E1 subset defined in scripts/e1_subsets.py: the full set as of 2026-11-20
(n=20, truncated tie group complete), the published 2026-11-17 subset (n=11)
and the active-only subset (n=6, pKd > 5.0).

The same null supplies an empirical 95th percentile of |rho| as a descriptive
random-ranking baseline. The p-value includes all permutations at least as
extreme as the observation, using exact integer rank dot products to avoid
floating-point loss of equal-tail events. Non-significance is not equivalence.

Reads only results/ext_ranking.csv (the point estimates for the published n=11
subset are re-derived and asserted against the numbers quoted in ext_results.md,
so this script fails loudly if the two ever drift apart). Writes
results/permutation.csv and results/permutation.md; deterministic under the fixed
seed (verified by running twice and comparing sha256).
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
sys.path.insert(0, str(HERE))
from e1_subsets import PUBLISHED_E1, subsets as e1_subsets  # noqa: E402

N_PERM = 20000
SEED = 42

# published point estimates (ext_results.md, E1 section) for the drift assert,
# which is checked on the 2026-11-17 subset so the number stays comparable
PUBLISHED = {"affinity_pred": -0.543, "vina_neg": -0.286}
ASSERT_SUBSET = "published_11"


def avgrank(v: np.ndarray) -> np.ndarray:
    return rankdata(v, method="average")


def perm_test(score: np.ndarray, pkd: np.ndarray, rng: np.random.Generator):
    """Two-sided Mantel permutation p for Spearman(score, pkd), pkd permuted."""
    n = len(pkd)
    rx = avgrank(score)
    ry = avgrank(pkd)
    rho_obs = float(spearmanr(score, pkd).statistic)
    keys = rng.random((N_PERM, n))
    perm_idx = np.argsort(keys, axis=1)
    Ry = ry[perm_idx]
    rx_z = (rx - rx.mean()) / rx.std()
    Ry_z = (Ry - Ry.mean(axis=1, keepdims=True)) / Ry.std(axis=1, keepdims=True)
    null = Ry_z @ rx_z / n
    abs_null = np.abs(null)
    # Average ranks are integers or half-integers. Twice-centered ranks are
    # integers, and their norm is invariant under permutation. Compare their
    # dot products exactly: separate floating-point correlation calculations
    # can otherwise exclude null values mathematically equal to the observation.
    rx2 = (2 * rx - (n + 1)).astype(np.int64)
    ry2 = (2 * ry - (n + 1)).astype(np.int64)
    obs_dot = abs(int(rx2 @ ry2))
    null_dot = np.abs(ry2[perm_idx] @ rx2)
    p_two = float((1 + (null_dot >= obs_dot).sum()) / (N_PERM + 1))
    return {"rho_obs": round(rho_obs, 4), "p_two_sided": p_two,
            "null_median_abs_rho": round(float(np.median(abs_null)), 4),
            "null_p95_abs_rho": round(float(np.percentile(abs_null, 95)), 4),
            "n_perm": N_PERM}


def nominal_association(score: str, rho: float, p: float) -> str:
    """Classify a two-sided, unadjusted association, not predictive utility."""
    expected_sign = {"affinity_pred": -1, "vina_neg": 1, "binary_bind_prob": 1}[score]
    if not (p < 0.05) or rho == 0:
        return "not_detected"
    return "expected_direction" if rho * expected_sign > 0 else "opposite_direction"


def main() -> None:
    e1 = pd.read_csv(RESULTS / "ext_ranking.csv")
    e1["vina_neg"] = -e1["vina_best"]
    scores = ["affinity_pred", "vina_neg", "binary_bind_prob"]
    groups = e1_subsets(e1)

    rows = []
    for key, (_label, sub) in groups.items():
        for col in scores:
            rng = np.random.default_rng(SEED)
            r = perm_test(sub[col].values.astype(float),
                          sub["true_pkd"].values.astype(float), rng)
            rows.append({"subset": key, "score": col, "n": len(sub), **r})

    tab = pd.DataFrame(rows)
    # drift asserts: the published-subset point estimates must match the report
    for col, pub in PUBLISHED.items():
        got = float(tab[(tab.subset == ASSERT_SUBSET) & (tab.score == col)].rho_obs.iloc[0])
        assert abs(got - pub) < 1e-3, f"{col}: recomputed {got} vs published {pub}"
    assert set(e1[e1["name"].isin(PUBLISHED_E1)]["name"]) == set(PUBLISHED_E1), \
        "the published 11 E1 runs are no longer all present in ext_ranking.csv"

    tab.to_csv(RESULTS / "permutation.csv", index=False, lineterminator="\n")

    lines = [
        "# B3 置换检验与随机基线（2026-11-18 首版，2026-11-20 随 E1 截断组补齐重跑，CPU，仅读已发布 ext_ranking.csv）",
        "",
        f"- 口径：Mantel 置换（整排 pKd 向量、保留 pKd=5.0 截断并列），双侧，{N_PERM} 次重排，"
        f"rng seed {SEED}，+1 校正；尾部计数用两倍中心化秩的整数内积，包含与观测值相等的置换，"
        "避免浮点误差漏计。与 ext_results.md 的 bootstrap CI 互补（同一数据、频率学派视角）。",
        "- 子集定义见 `scripts/e1_subsets.py`：`all_e1` 为当前全部 E1 配体，"
        f"`published_11` 为 2026-11-17 首批（6 活性 + Vina 最强 5 截断，与既有 README 引用同源），"
        "`active_pkd_gt_5` 为仅活性子集。",
    ]
    for key, (label, sub) in groups.items():
        lines += ["", f"## {label}（n={len(sub)}）", "",
                  "| 分数 | Spearman ρ | 置换 p（双侧） | 零分布\\|ρ\\|中位 | 零分布\\|ρ\\| 95分位 |",
                  "|---|---|---|---|---|"]
        for col in scores:
            r = tab[(tab.subset == key) & (tab.score == col)].iloc[0]
            lines.append(f"| {col} | {r.rho_obs:+.3f} | {r.p_two_sided:.4f} | "
                         f"{r.null_median_abs_rho:.3f} | {r.null_p95_abs_rho:.3f} |")
        for col in ("affinity_pred", "vina_neg"):
            r = tab[(tab.subset == key) & (tab.score == col)].iloc[0]
            state = nominal_association(col, r.rho_obs, r.p_two_sided)
            reading = {
                "expected_direction": "未校正、方向正确的关联信号",
                "opposite_direction": "未校正、方向相反的关联信号",
                "not_detected": "未检出关联信号（不等于证明随机或等效）",
            }[state]
            lines.append(f"- {col}：ρ={r.rho_obs:+.3f}，未校正 p={r.p_two_sided:.4f}；"
                         f"零分布 |ρ| 的经验 95 分位为 {r.null_p95_abs_rho:.3f}，"
                         f"观测 |ρ|={abs(r.rho_obs):.3f} ⇒ {reading}。")
    aff11 = tab[(tab.subset == "published_11") & (tab.score == "affinity_pred")].iloc[0]
    lines += [
        "",
        "## 分辨率与既有结论的关系",
        "",
        f"- n=6 有 6!=720 种索引排列；本脚本仍使用 {N_PERM} 次 Monte Carlo 重排和 +1 校正，"
        f"数值步长为 1/({N_PERM}+1)，不是精确枚举的 1/720，也不是 1/721。"
        "重复排列、并列和双侧尾部会影响实际可达 p 值；增加重排次数不会增加配体数。",
        "- affinity_pred 是 log10(IC50/μM)，低值对应强结合，与 pKd 的预期相关为负；"
        "vina_neg 与 binary_bind_prob 的预期方向为正。",
        f"- n=11 的 Boltz-2 ρ={aff11.rho_obs:+.3f} 与 2026-11-18 首版一致（脚本保留 ρ 的漂移断言）。"
        f"2026-11-20 修复等值置换漏计后，p 由 0.0848 更正为 {aff11.p_two_sided:.4f}；"
        "原始预测、配体集合、相关系数和 bootstrap CI 未改变，p 值更正不属于模型行为变化。",
        "- 补齐后的全量集合只是把截断并列组从『被挑过的一半』换成整组；结论方向若发生变化，"
        "必须同时说明这是截断组构成变化所致，而不是模型行为变化。",
        f"- 本表共 {len(tab)} 个分数×子集比较，p 值均未做多重比较校正。即便某个子集出现"
        "方向正确的名义关联，也不能据此证明预测可用性，或宣称 Boltz-2 优于 Vina。"
        "全量的活性/截断组分离也不能替代仅活性子集的精细排序能力。",
        "",
        "复跑：`python scripts/permutation_test.py`（仅 pandas/numpy/scipy；连跑两次产物哈希一致）。",
    ]
    (RESULTS / "permutation.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print("\n".join(lines))


if __name__ == "__main__":
    main()
