"""B3: permutation tests and random-ranking baselines for the E1 ranking result.

The E1 verdict in ext_results.md rests on bootstrap CIs spanning zero. This
script adds the complementary frequentist view. Under the null "the score
carries no information about pKd", Spearman rho is exactly distributed over
the pairings of the two rank vectors (Mantel permutation test). Ligands share
the DAVIS truncation value pKd = 5.0, so the pKd vector is permuted as a whole
against the fixed score ranks, which preserves the tie structure.

Monte-Carlo two-sided p-values (20000 shuffles, rng seed 42, +1 correction) for
the three E1 scores (Boltz-2 affinity, -Vina, binary bind probability) on every
E1 subset defined in scripts/e1_subsets.py: the full set as of 2026-11-20
(n=20, truncated tie group complete), the published 2026-11-17 subset (n=11)
and the active-only subset (n=6, pKd > 5.0).

The same null doubles as the random baseline: the |rho| a coin-flip ranking
reaches 95% of the time. If the observed |rho| sits below that percentile,
the data cannot distinguish the model from random ordering even before
considering the sign.

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
    p_two = float((1 + (abs_null >= abs(rho_obs)).sum()) / (N_PERM + 1))
    return {"rho_obs": round(rho_obs, 4), "p_two_sided": p_two,
            "null_median_abs_rho": round(float(np.median(abs_null)), 4),
            "null_p95_abs_rho": round(float(np.percentile(abs_null, 95)), 4),
            "n_perm": N_PERM}


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
        f"rng seed {SEED}，+1 校正；与 ext_results.md 的 bootstrap CI 互补（同一数据、频率学派视角）。",
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
            rnd = "与随机排序不可区分" if abs(r.rho_obs) < r.null_p95_abs_rho else "超出随机 95 分位"
            lines.append(f"- {col}：ρ={r.rho_obs:+.3f}，p={r.p_two_sided:.4f}"
                         f"{'（<0.05）' if r.p_two_sided < 0.05 else '（0.05 水平不显著）'}；"
                         f"随机排序在 n={len(sub)} 下有 5% 概率碰出 |ρ|≥{r.null_p95_abs_rho:.3f}，"
                         f"观测 |ρ|={abs(r.rho_obs):.3f} ⇒ {rnd}。")
    aff11 = tab[(tab.subset == "published_11") & (tab.score == "affinity_pred")].iloc[0]
    lines += [
        "",
        "## 分辨率与既有结论的关系",
        "",
        "- n=6 时只有 6!=720 种配对，置换 p 的分辨率下限约 1/721≈0.0014；该子集只能给出"
        "『显著/不显著』的粗判，不宜作精细排序结论。",
        f"- n=11 的 Boltz-2 结果（ρ={aff11.rho_obs:+.3f}，p={aff11.p_two_sided:.4f}）与 2026-11-18 首版一致"
        "（脚本对该子集的 ρ 做了漂移断言，数值变了就会报错而不是悄悄改口径）。",
        "- 补齐后的全量集合只是把截断并列组从『被挑过的一半』换成整组；结论方向若发生变化，"
        "必须同时说明这是截断组构成变化所致，而不是模型行为变化。",
        "- 不得写成『Boltz-2 显著差于随机』或『显著优于随机』——两个方向都没有证据。",
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
