"""Sequence extraction from the crystal PDB (vendored from the local probe).

This is the original probe/make_yaml.py with ONLY the path constants changed to
repository-relative ones (see README, 复现 section). The logic — construct
sequence extraction, the CHOICE table of six baseline pairs, and the emitted
YAML (including the `msa: empty` line that later turned out to make the
baseline runs single-sequence) — is byte-identical to what produced the
published results.

Original docstring kept verbatim below.

Boltz-2 feasibility probe: six protein-ligand pairs taken from the rows that
AutoDock Vina already scored in dti-affinity-bench, so the three-way table
(true pKd / Vina / Boltz-2) is on identical rows and nothing has to be re-docked.

Two deliberate design choices:

1. The protein input is the *crystal construct* sequence extracted from the very
   PDB the docking arm used (8c7x / 7ukz / 6yz4 chain A), not the full-length
   DAVIS UniProt sequence. That keeps the comparison honest (same physical
   object Vina was handed) and keeps token counts at 274-296 instead of 766,
   which is itself part of what we are measuring on a 10.6 GB free GPU.
2. Ligands are picked to span measured affinity *within one target*, so the
   affinity head is graded on ranking, not on absolute values:
   BRAF 9.73 vs 6.43, CDC2L1 6.54 vs 5.89, MKK7 7.52 vs 5.28.

Everything is read from the benchmark repo; no SMILES or sequence is retyped.
"""

from __future__ import annotations

import os
import pathlib

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCH = pathlib.Path(os.environ.get("DTI_BENCH", str(ROOT / "data" / "dti_subset")))
PANEL = BENCH / "results" / "docking" / "target_panel.csv"
PAIRS_CSV = BENCH / "results" / "docking" / "pairs.csv"
LIGANDS = BENCH / "data" / "processed" / "ligands.csv"
STRUCT_DIR = BENCH / "data" / "raw" / "docking"
OUTDIR = ROOT / "inputs"

# target_id -> ligand_ids, all present in results/docking/pairs.csv
CHOICE = {
    "T0013": ["L02069", "L00585"],   # BRAF   9.73 / 6.43
    "T0004": ["L01091", "L01149"],   # CDC2L1 6.54 / 5.89
    "T0015": ["L00602", "L00770"],   # MKK7   7.52 / 5.28
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # modified / modelled residues that appear in these three entries
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "ALY": "K",
    "MLZ": "K", "NLE": "X", "AIB": "X",
}


def construct_sequence(pdb_path: pathlib.Path, chain: str) -> tuple[str, int, int]:
    """One-letter sequence of `chain` from ATOM records, plus gaps and drops.

    PDB column layout: residue name 18-20, chain 22, residue number 23-26
    (0-based slice [17:20], [21], [22:26]).  HETATM records are ignored on
    purpose -- waters, EDO cryoprotectant and the co-crystal ligand must not
    end up in the polymer sequence.
    """
    seen: dict[int, str] = {}
    dropped: list[str] = []
    with pdb_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue
            if line[21] != chain:
                continue
            name = line[17:20].strip().upper()
            one = THREE_TO_ONE.get(name)
            if one is None:
                dropped.append(f"{name}{line[22:26].strip()}")
                continue
            seen.setdefault(int(line[22:26]), one)
    numbers = sorted(seen)
    gaps = sum(1 for a, b in zip(numbers, numbers[1:]) if b - a != 1)
    return "".join(seen[n] for n in numbers), gaps, len(dropped)


def main() -> None:
    panel = pd.read_csv(PANEL).set_index("target_id")
    docked = pd.read_csv(PAIRS_CSV)
    lig = pd.read_csv(LIGANDS).set_index("ligand_id")

    OUTDIR.mkdir(exist_ok=True)
    meta: list[dict[str, object]] = []
    for tid, lids in CHOICE.items():
        pdb_id = str(panel.at[tid, "pdb_id"]).lower()
        chain = str(panel.at[tid, "chain"])
        pdb_path = STRUCT_DIR / f"{pdb_id}.pdb"
        seq, gaps, dropped = construct_sequence(pdb_path, chain)
        gene = str(panel.at[tid, "gene"])
        for lid in lids:
            rows = docked[(docked.target_id == tid) & (docked.ligand_id == lid)]
            if rows.empty:
                raise SystemExit(f"{tid}/{lid} is not in the docked set")
            smi = str(lig.at[lid, "canonical_smiles"])
            if "." in smi:  # multi-component salts: keep the largest fragment
                smi = max(smi.split("."), key=len)
            name = f"{tid}_{lid}"
            (OUTDIR / f"{name}.yaml").write_text(
                "version: 1\n"
                "sequences:\n"
                "  - protein:\n"
                f"      id: A\n"
                f"      sequence: {seq}\n"
                # boltz refuses a protein chain with neither an MSA file nor
                # --use_msa_server; "empty" is the documented single-sequence mode
                "      msa: empty\n"
                "  - ligand:\n"
                "      id: B\n"
                f"      smiles: {smi}\n"
                "properties:\n"
                "  - affinity:\n"
                "      binder: B\n",
                encoding="utf-8",
            )
            meta.append(
                {
                    "name": name,
                    "target_id": tid,
                    "gene": gene,
                    "pdb_id": panel.at[tid, "pdb_id"],
                    "construct_len": len(seq),
                    "chain_breaks": gaps,
                    "dropped_residues": dropped,
                    "ligand_id": lid,
                    "n_docked_rows": len(rows),
                    "true_pkd": round(float(rows.pkd.mean()), 3),
                    "vina_best": round(float(rows.vina_best.mean()), 3),
                    "smiles": smi,
                }
            )
    df = pd.DataFrame(meta)
    (ROOT / "results" / "probe_meta.csv").parent.mkdir(exist_ok=True)
    df.to_csv(ROOT / "results" / "probe_meta.csv", index=False)
    print(f"wrote {len(df)} systems to {OUTDIR}")
    print(df[["name", "gene", "construct_len", "chain_breaks", "dropped_residues",
              "true_pkd", "vina_best"]].to_string(index=False))


if __name__ == "__main__":
    main()
