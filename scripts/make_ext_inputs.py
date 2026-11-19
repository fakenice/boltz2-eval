"""Generate inputs for the Boltz-2 extension experiments (2026-11-17).

E1  BRAF (T0013) within-target ranking expansion: every ligand with measured
    pKd > 5.0 (6 rows) plus the 5 pKd == 5.0 rows with the strongest Vina
    scores (adversarial: docking thinks they are tight, measurement says weak).
    All runs reuse the persisted BRAF a3m so the MSA is identical across ligands.
E2  MSA ablation, corrected protocol (2026-11-17 finding): the original six
    probe/out/msa_* runs carried `msa: empty` in their YAML, which parses to
    msa_id=-1 and takes precedence over --use_msa_server (schema.py L1118,
    main.py L517-529); their records confirm msa_id=-1. So the originals are
    the SINGLE-SEQUENCE baseline despite the "+msa" preset label. The ablation
    therefore adds real-MSA runs (persisted local a3m) for the four non-BRAF
    originals; the two BRAF real-MSA points come from E1. One single-sequence
    re-run (T0013_L00585) checks that the original numbers reproduce.
E3  Sampling-quality curve at working size: two complexes at two structure
    presets (25 steps/1 recycle and 100 steps/2 recycles), affinity head kept
    at its default (5 x 200) in every run so only structure cost varies. The
    200/3 point of each curve reuses the real-MSA default run (E1_T0013_L02069
    / E2_T0004_L01091) so the whole curve is one MSA protocol.
EQ  Protocol equivalence: one run with msa auto + --use_msa_server (the path
    the original report believed it had used) vs the local-a3m E1 run of the
    same complex, to verify persisted-MSA reuse matches the server path.

Nothing is retyped: sequences come from the crystal PDB chain A via
make_yaml.construct_sequence (vendored in scripts/), SMILES and truth from the
bundled benchmark subset under data/dti_subset.

Compared with the local working copy, only path constants changed (repo-
relative, with DTI_BENCH env override); the E1 docstring count was corrected
from "11 rows" to the actual 6 actives + 5 truncated = 11 runs. The published
inputs/*.yaml carry the a3m as a repo-root-relative path ("msa/T0013.a3m",
run from the repo root); the original runs used absolute paths to the byte-
identical files under msa/, which are published here unchanged.
"""
from __future__ import annotations

import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent          # scripts/
ROOT = HERE.parent                                       # repo root
sys.path.insert(0, str(HERE))
from make_yaml import construct_sequence, THREE_TO_ONE  # noqa: E402,F401

import pandas as pd  # noqa: E402

BENCH = pathlib.Path(os.environ.get("DTI_BENCH", str(ROOT / "data" / "dti_subset")))
PANEL = BENCH / "results" / "docking" / "target_panel.csv"
PAIRS = BENCH / "results" / "docking" / "pairs.csv"
LIGANDS = BENCH / "data" / "processed" / "ligands.csv"
STRUCT_DIR = BENCH / "data" / "raw" / "docking"
INDIR = ROOT / "inputs"

ORIGINAL = {
    "T0013": ["L02069", "L00585"],
    "T0004": ["L01091", "L01149"],
    "T0015": ["L00602", "L00770"],
}


def yaml_text(seq: str, smi: str, msa: str) -> str:
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {seq}\n"
        f"      msa: {msa}\n"
        "  - ligand:\n"
        "      id: B\n"
        f"      smiles: {smi}\n"
        "properties:\n"
        "  - affinity:\n"
        "      binder: B\n"
    )


def main() -> None:
    panel = pd.read_csv(PANEL).set_index("target_id")
    pairs = pd.read_csv(PAIRS)
    lig = pd.read_csv(LIGANDS).set_index("ligand_id")
    INDIR.mkdir(exist_ok=True)

    seqs = {}
    for tid in ORIGINAL:
        pdb_id = str(panel.at[tid, "pdb_id"]).lower()
        chain = str(panel.at[tid, "chain"])
        seq, gaps, dropped = construct_sequence(STRUCT_DIR / f"{pdb_id}.pdb", chain)
        seqs[tid] = seq
        print(f"{tid} {pdb_id} chain {chain}: {len(seq)} aa, gaps={gaps}, dropped={dropped}")

    def smiles(lid: str) -> str:
        smi = str(lig.at[lid, "canonical_smiles"])
        return max(smi.split("."), key=len) if "." in smi else smi

    meta = []

    def add(name, exp, tid, lid, msa, msa_mode=None):
        rows = pairs[(pairs.target_id == tid) & (pairs.ligand_id == lid)]
        assert len(rows) >= 1, f"{tid}/{lid} missing from docked set"
        (INDIR / f"{name}.yaml").write_text(yaml_text(seqs[tid], smiles(lid), msa), encoding="utf-8")
        if msa_mode is None:
            msa_mode = "empty" if msa == "empty" else "local_a3m"
        meta.append({
            "name": name, "exp": exp, "target_id": tid,
            "gene": panel.at[tid, "gene"], "pdb_id": panel.at[tid, "pdb_id"],
            "construct_len": len(seqs[tid]), "ligand_id": lid,
            "true_pkd": round(float(rows.pkd.mean()), 3),
            "vina_best": round(float(rows.vina_best.mean()), 3),
            "msa_mode": msa_mode,
        })

    # E1: BRAF ranking expansion (pairs.csv can hold >1 row per ligand across
    # published parts; aggregate per ligand first, matching add()'s mean semantics)
    braf = (pairs[pairs.target_id == "T0013"]
            .groupby("ligand_id", as_index=False)
            .agg(pkd=("pkd", "mean"), vina_best=("vina_best", "mean")))
    active = braf[braf.pkd > 5.0]
    trunc_all = braf[braf.pkd == 5.0].sort_values("vina_best")
    # The 2026-11-17 batch ran only the 5 strongest-Vina truncated rows (the
    # adversarial subset). trunc_rest are the remaining 9 pKd == 5.0 ligands,
    # added 2026-11-20 so the truncated group is complete rather than curated;
    # they stay in the same E1 group because the question is "can Boltz-2 rank
    # within BRAF at all", and a hand-picked half of a tie group is an answer
    # that cannot be re-asked. Analyses must keep reporting the published
    # n=11 subset separately (see scripts/permutation_test.py).
    trunc = trunc_all.head(5)
    trunc_rest = trunc_all.iloc[5:]
    sel = pd.concat([active, trunc, trunc_rest])
    print(f"E1 selection: {len(active)} active + {len(trunc)} truncated (2026-11-17) "
          f"+ {len(trunc_rest)} truncated (2026-11-20) = {len(sel)}")
    for r in sel.itertuples():
        add(f"E1_T0013_{r.ligand_id}", "E1", "T0013", r.ligand_id, "msa/T0013.a3m")

    # E2: real-MSA side of the ablation for the four non-BRAF originals (the
    # two BRAF real-MSA points are E1_T0013_L02069 / E1_T0013_L00585), plus one
    # single-sequence re-run to check the original numbers reproduce.
    for tid, lids in [("T0004", ["L01091", "L01149"]), ("T0015", ["L00602", "L00770"])]:
        for lid in lids:
            add(f"E2_{tid}_{lid}", "E2", tid, lid, f"msa/{tid}.a3m")
    add("E2rep_T0013_L00585", "E2rep", "T0013", "L00585", "empty")

    # E3: sampling curve (two complexes, two structure presets; affinity default)
    for tid, lid in [("T0013", "L02069"), ("T0004", "L01091")]:
        for tag in ("E3a", "E3b"):
            add(f"{tag}_{tid}_{lid}", tag, tid, lid, f"msa/{tid}.a3m")

    # EQ: server-MSA path (msa auto + --use_msa_server) vs local-a3m reuse
    add("EQ_T0013_L02069", "EQ", "T0013", "L02069", "0", msa_mode="server")

    df = pd.DataFrame(meta)
    (ROOT / "results").mkdir(exist_ok=True)
    df.to_csv(ROOT / "results" / "ext_meta.csv", index=False)
    print(f"wrote {len(df)} inputs")
    print(df.groupby("exp").size().to_string())


if __name__ == "__main__":
    main()
