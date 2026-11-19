"""Shared E1 subset definitions (added 2026-11-20).

E1 asks one question: can Boltz-2 rank ligands within BRAF? The 2026-11-17 batch
answered it on 11 ligands — every measured active (pKd > 5.0) plus the 5 truncated
rows (pKd == 5.0, the DAVIS Kd ceiling) with the strongest Vina scores, i.e. the
subset docking is most confident about. The 2026-11-20 batch added the remaining 9
truncated ligands so the tie group is complete instead of curated.

Both views have to stay reportable: the published n=11 numbers are quoted in
README.md, results/ext_results.md and results/permutation.md, while the full set is
the one that answers the question without a hand-picked half of a tie group. The
subsets are therefore derived from the data by the documented selection rule *and*
pinned by name here, so a change in the rule or in the run set fails loudly instead
of silently re-labelling a published number.

Used by scripts/analyze_ext.py (E1 ranking block) and scripts/permutation_test.py.
"""
from __future__ import annotations

import pandas as pd

# The 11 runs of the 2026-11-17 batch, by ext_runs.csv `name`.
PUBLISHED_E1 = (
    "E1_T0013_L02069", "E1_T0013_L00585", "E1_T0013_L00349", "E1_T0013_L00772",
    "E1_T0013_L00808", "E1_T0013_L01163", "E1_T0013_L00602", "E1_T0013_L01149",
    "E1_T0013_L00261", "E1_T0013_L01299", "E1_T0013_L01320",
)
# pKd == 5.0 is the DAVIS truncation ceiling; above it a ligand is measured active.
TRUNCATION = 5.0


def actives(e1: pd.DataFrame) -> pd.DataFrame:
    return e1[e1.true_pkd > TRUNCATION]


def truncated(e1: pd.DataFrame) -> pd.DataFrame:
    return e1[e1.true_pkd == TRUNCATION]


def published_11(e1: pd.DataFrame) -> pd.DataFrame:
    """Re-derive the 2026-11-17 subset from whatever E1 rows are present.

    With only the published rows on disk this is the identity; with the full
    truncated tie group it must still select the same 11 names (the 9 additions are
    strictly weaker by Vina), which is what the assert pins. Rows are selected by
    name mask so the input's row order is preserved: the Monte-Carlo permutation
    draws depend on that order, and re-ordering would quietly change the published
    p-values while leaving every observed rho alone.
    """
    keep = set(actives(e1)["name"]) | set(truncated(e1).sort_values("vina_best").head(5)["name"])
    sub = e1[e1["name"].isin(keep)]
    names = set(sub["name"])
    assert names == set(PUBLISHED_E1), (
        f"E1 published subset re-derived as {sorted(names)}; published batch is "
        f"{sorted(PUBLISHED_E1)}. The selection rule or the run set changed — "
        "update the docs that quote n=11 before trusting this.")
    return sub


def subsets(e1: pd.DataFrame) -> dict[str, tuple[str, pd.DataFrame]]:
    """Ordered {key: (human label, rows)} for the E1 ranking analyses."""
    out = {
        "all_e1": ("全部 E1 配体", e1),
        "published_11": ("2026-11-17 已发布子集（6 活性 + Vina 最强 5 截断）", published_11(e1)),
        "active_pkd_gt_5": (f"仅活性子集（pKd>{TRUNCATION:g}）", actives(e1)),
    }
    return {k: v for k, v in out.items() if len(v[1]) >= 4}
