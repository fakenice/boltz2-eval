"""B3: permutation tests and random-ranking baselines for the E1 ranking result.

The E1 verdict in ext_results.md rests on bootstrap CIs spanning zero. This
script adds the complementary frequentist view. Under the null "the score
carries no information about pKd", Spearman rho is exactly distributed over
the pairings of the two rank vectors (Mantel permutation test). Five ligands
share the DAVIS truncation value pKd = 5.0, so the pKd vector is permuted as a
whole against the fixed score ranks, which preserves the tie structure.

Monte-Carlo two-sided p-values (20000 shuffles, rng seed 42, +1 correction)
for the three E1 scores (Boltz-2 affinity, -Vina, binary bind probability) on
both the full set (n=11) and the active-only subset (n=6, pKd > 5.0).

The same null doubles as the random baseline: the |rho| a coin-flip ranking
reaches 95% of the time. If the observed |rho| sits below that percentile,
the data cannot distinguish the model from random ordering even before
considering the sign.

Reads only results/ext_ranking.csv (published; the point estimates are
re-derived and asserted against the numbers quoted in ext_results.md, so
this script fails loudly if the two ever drift apart). Writes
results/permutation.csv and results/permutation.md; deterministic under the
fixed seed (verified by running twice and comparing sha256).
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"

N_PERM = 20000
SEED = 42

# published point estimates (ext_results.md, E1 section) for the drift assert
PUBLISHED = {"affinity_pred": -0.543, "vina_neg": -0.286}


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
    subsets = {"all_11": e1, "active_pkd_gt_5": e1[e1.true_pkd > 5.0]}
    scores = ["affinity_pred", "vina_neg", "binary_bind_prob"]

    rows = []
    for sname, sub in subsets.items():
        for col in scores:
            rng = np.random.default_rng(SEED)
            r = perm_test(sub[col].values.astype(float),
                          sub["true_pkd"].values.astype(float), rng)
            rows.append({"subset": sname, "score": col, "n": len(sub), **r})

    tab = pd.DataFrame(rows)
    # drift asserts: the n=11 point estimates must match the published E1 numbers
    for col, pub in PUBLISHED.items():
        got = float(tab[(tab.subset == "all_11") & (tab.score == col)].rho_obs.iloc[0])
        assert abs(got - pub) < 1e-3, f"{col}: recomputed {got} vs published {pub}"

    tab.to_csv(RESULTS / "permutation.csv", index=False, lineterminator="\n")

    aff = tab[(tab.subset == "all_11") & (tab.score == "affinity_pred")].iloc[0]
    vin = tab[(tab.subset == "all_11") & (tab.score == "vina_neg")].iloc[0]
    lines = [
        "# B3 置换检验与随机基线（2026-11-18，CPU，仅读已发布 ext_ranking.csv）",
        "",
        f"- 口径：Mantel 置换（整排 pKd 向量、保留 5 个 pKd=5.0 截断并列），双侧，{N_PERM} 次重排，"
        f"rng seed {SEED}，+1 校正；与 ext_results.md 的 bootstrap CI 互补（同一数据、频率学派视角）。",
        "",
        "## 全部 11 配体（n=11）",
        "",
        "| 分数 | Spearman ρ | 置换 p（双侧） | 零分布\\|ρ\\|中位 | 零分布\\|ρ\\| 95分位 |",
        "|---|---|---|---|---|",
    ]
    for col in scores:
        r = tab[(tab.subset == "all_11") & (tab.score == col)].iloc[0]
        lines.append(f"| {col} | {r.rho_obs:+.3f} | {r.p_two_sided:.4f} | "
                     f"{r.null_median_abs_rho:.3f} | {r.null_p95_abs_rho:.3f} |")
    lines += [
        "",
        f"- Boltz-2 亲和力：ρ={aff.rho_obs:+.3f}，p={aff.p_two_sided:.4f}"
        f"{'（<0.05，但方向为负：排序系统性反了，不是「没有信号」）' if aff.p_two_sided < 0.05 else '（0.05 水平不显著）'}；"
        f"随机排序在 n=11 下也有 5% 概率碰出 |ρ|≥{aff.null_p95_abs_rho:.3f}，观测 |ρ|={abs(aff.rho_obs):.3f} "
        f"{'低于' if abs(aff.rho_obs) < aff.null_p95_abs_rho else '高于'}该线 ⇒ "
        f"{'与随机排序不可区分' if abs(aff.rho_obs) < aff.null_p95_abs_rho else '超出随机 95 分位'}。",
        f"- Vina：ρ={vin.rho_obs:+.3f}，p={vin.p_two_sided:.4f}"
        f"{'，观测 |ρ| 低于随机 95 分位 ⇒ 与随机不可区分' if abs(vin.rho_obs) < vin.null_p95_abs_rho else '，观测 |ρ| 高于随机 95 分位'}。",
        "",
        "## 仅活性子集（pKd>5.0，n=6）",
        "",
        "| 分数 | Spearman ρ | 置换 p（双侧） | 零分布\\|ρ\\|中位 | 零分布\\|ρ\\| 95分位 |",
        "|---|---|---|---|---|",
    ]
    for col in scores:
        r = tab[(tab.subset == "active_pkd_gt_5") & (tab.score == col)].iloc[0]
        lines.append(f"| {col} | {r.rho_obs:+.3f} | {r.p_two_sided:.4f} | "
                     f"{r.null_median_abs_rho:.3f} | {r.null_p95_abs_rho:.3f} |")
    lines += [
        "",
        "n=6 时 6!=720 种配对全部枚举也才 720，置换 p 的分辨率下限约 1/721≈0.0014；"
        "该子集只能给出『显著/不显著』的粗判，不宜作精细排序结论。",
        "",
        "## 与 E1 bootstrap 结论的关系",
        "",
        "- E1 报 CI 跨 0（排序能力未建立）；本脚本给出同一事实的 p 值表述，并补上『随机基线』：",
        f"  n=11 时随机排序的 |ρ| 95 分位：Boltz-2 零分布为 {aff.null_p95_abs_rho:.3f}"
        f"（观测 |ρ|={abs(aff.rho_obs):.3f}），Vina 零分布为 {vin.null_p95_abs_rho:.3f}"
        f"（观测 |ρ|={abs(vin.rho_obs):.3f}）。",
        "- 不得写成『Boltz-2 显著差于随机』或『显著优于随机』——两个方向都没有证据。",
        "",
        "复跑：`python scripts/permutation_test.py`（仅 pandas/numpy/scipy；连跑两次产物哈希一致）。",
    ]
    (RESULTS / "permutation.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    try:
        import sys
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print("\n".join(lines))


if __name__ == "__main__":
    main()
