"""Fetch MMSeqs2 MSA once per construct and save a3m for reuse across runs.

Same protocol as boltz --use_msa_server (boltz.data.msa.mmseqs2.run_mmseqs2,
mode=env, filter on, no pairing), but the a3m is persisted so every ligand run
on the same construct consumes an identical MSA and needs no network.
Sequences are re-extracted from the crystal PDB chain A exactly as make_yaml.py
did for the original six complexes.

The published msa/*.a3m were produced by this script; re-running only matters
if you want fresh retrievals (note ColabFold server results can evolve).
Requires the boltz package (pip install boltz==2.0.3) and network access.
Only path constants changed vs the local working copy.
"""
from __future__ import annotations

import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from make_yaml import construct_sequence  # noqa: E402

import pandas as pd  # noqa: E402
from boltz.data.msa.mmseqs2 import run_mmseqs2  # noqa: E402

BENCH = pathlib.Path(os.environ.get("DTI_BENCH", str(ROOT / "data" / "dti_subset")))
PANEL = BENCH / "results" / "docking" / "target_panel.csv"
STRUCT_DIR = BENCH / "data" / "raw" / "docking"

TARGETS = ["T0013", "T0004", "T0015"]  # BRAF 8c7x / CDC2L1 7ukz / MKK7 6yz4


def main() -> None:
    panel = pd.read_csv(PANEL).set_index("target_id")
    (ROOT / "msa").mkdir(exist_ok=True)
    for tid in TARGETS:
        pdb_id = str(panel.at[tid, "pdb_id"]).lower()
        chain = str(panel.at[tid, "chain"])
        seq, gaps, dropped = construct_sequence(STRUCT_DIR / f"{pdb_id}.pdb", chain)
        out = ROOT / "msa" / f"{tid}.a3m"
        if out.exists() and out.stat().st_size > 0:
            print(f"{tid}: cached {out.stat().st_size} bytes", flush=True)
            continue
        print(f"{tid} ({pdb_id} chain {chain}, {len(seq)} aa): submitting", flush=True)
        a3m_lines = run_mmseqs2(seq, prefix=str(ROOT / "msa" / f"{tid}_work"))
        out.write_text(a3m_lines[0], encoding="utf-8")
        n_seqs = sum(1 for ln in a3m_lines[0].splitlines() if ln.startswith(">"))
        print(f"{tid}: wrote {out.name}, {n_seqs} MSA sequences, {out.stat().st_size} bytes", flush=True)


if __name__ == "__main__":
    main()
