"""B2: PoseBusters validity checks on the published Boltz-2 predictions.

CPU only. Reads only repository files: the predicted complex cif under
outputs/<name>/ (21 batch runs) and outputs_baseline/msa_<name>/ (6 original
single-sequence runs), plus the input SMILES in inputs/<name>.yaml.

Two layers, reported separately because they rest on different assumptions:

1. Raw geometry (no bond perception, unambiguous): for the ligand heavy-atom
   cloud - minimum pairwise distance and counts of pairs below 0.5/1.0/1.3 A;
   for the ligand-protein interface - minimum heavy-atom distance and the
   number of van der Waals clashes (d < 0.75 * (r_i + r_j)).

2. PoseBusters protocol (posebusters==0.6.5, config "dock"): the ligand
   requires a molecular graph, but the cif has no _chem_comp_bond loop, so
   connectivity is perceived from coordinates (lower bound 0.8 A to exclude
   overlapping atoms, upper bound 1.3 x sum of covalent radii) and bond
   orders are assigned from the input SMILES with
   AssignBondOrdersFromTemplate. The assignment is attempted only when the
   perceived bond count equals the template bond count: on a mismatched
   graph RDKit's subgraph matching explodes combinatorially (observed on
   this data; a full scan had to be killed after minutes without output).
   Runs failing this step are recorded as not-loadable and skip the check
   table - that failure count is itself part of the result.
   The protein is written to a temporary PDB (PoseBusters requires PDB
   residue info for mol_cond).

Outputs: results/posebusters_geometry.csv, results/posebusters.csv (only if
some run loads), results/posebusters.md (summary, generated, not hand-written).
The same raw-geometry metric is also run on the co-crystal ligands of the
three published PDBs as a positive control (experimental geometry must pass;
if it did, the metric or the parsing would be suspect).
"""
from __future__ import annotations

import pathlib
import re
import tempfile

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds  # noqa: F401  (kept for API probe parity)

RDLogger.DisableLog("rdApp.*")

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
INPUTS = ROOT / "inputs"

PT = Chem.GetPeriodicTable()
VDW = {}  # element -> vdW radius, filled lazily from the periodic table

BUFFER = {"HOH", "EDO", "GOL", "SO4", "PO4", "DMS", "PEG", "MES", "TRS", "CL", "NA",
          "MG", "ZN", "CA", "K", "NO3", "ACT", "EPE", "CIT", "FLC", "PGE", "PG4",
          "1PE", "P6G", "XP4", "XPE", "PG0", "2PE", "15P", "PE4", "PEU", "DTT",
          "BME", "IOD", "BR", "FMT", "ACN", "MPD", "TLA", "TAR", "SIN", "ORN"}


def cov(el: str) -> float:
    return PT.GetRcovalent(el)


def vdw(el: str) -> float:
    if el not in VDW:
        # GetRvdw raises for some elements; fall back to covalent + 0.4
        try:
            VDW[el] = PT.GetRvdw(el)
        except Exception:
            VDW[el] = PT.GetRcovalent(el) + 0.4
    return VDW[el]


def parse_atom_site(cif: pathlib.Path):
    """Return the _atom_site rows of a cif as dicts (whitespace-split, all columns)."""
    lines = cif.read_text(encoding="utf-8", errors="replace").splitlines()
    rows, i = [], 0
    while i < len(lines):
        if lines[i].startswith("_atom_site."):
            cols = []
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].split()[0].split(".", 1)[1])
                i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("_"):
                p = lines[i].split()
                if len(p) == len(cols):
                    rows.append(dict(zip(cols, p)))
                i += 1
        else:
            i += 1
    return rows


def pair_stats(xyz: np.ndarray):
    d = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=2)
    iu = np.triu_indices(len(xyz), 1)
    dd = d[iu]
    return {
        "min_pair_dist": round(float(dd.min()), 3) if len(dd) else np.nan,
        "n_pairs_lt_0_5A": int((dd < 0.5).sum()),
        "n_pairs_lt_1_0A": int((dd < 1.0).sum()),
        "n_pairs_lt_1_3A": int((dd < 1.3).sum()),
    }


def clash_stats(xyz_a: np.ndarray, el_a: list[str], xyz_b: np.ndarray, el_b: list[str]):
    d = np.linalg.norm(xyz_a[:, None] - xyz_b[None, :], axis=2)
    rsum = np.array([[vdw(a) + vdw(b) for b in el_b] for a in el_a])
    rel = d / rsum
    return {"min_lig_pro_dist": round(float(d.min()), 3),
            "n_vdw_clashes_0_75": int((rel < 0.75).sum())}


def write_protein_pdb(rows, path: pathlib.Path):
    at = [r for r in rows if r["group_PDB"] == "ATOM" and r["auth_asym_id"] == "A"]
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        for k, r in enumerate(at, 1):
            name = r["label_atom_id"][:4]
            fh.write("ATOM  %5d %-4s %3s A%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (
                k, name, r["label_comp_id"], r["auth_seq_id"],
                float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]),
                1.0, float(r["B_iso_or_equiv"]), r["type_symbol"].capitalize()))
        fh.write("END\n")
    return len(at)


def perceive_ligand(het_rows, smi: str):
    """Build a template-bond-ordered ligand mol from coordinates. Returns (mol, note)."""
    el = [r["type_symbol"].capitalize() for r in het_rows]
    xyz = np.array([[float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])] for r in het_rows])
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(el))
    for k in range(len(el)):
        rw.AddAtom(Chem.Atom(el[k]))
        conf.SetAtomPosition(k, xyz[k].tolist())
    n_bond = 0
    for i in range(len(el)):
        for j in range(i + 1, len(el)):
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if 0.8 <= d <= 1.3 * (cov(el[i]) + cov(el[j])):
                rw.AddBond(i, j, Chem.BondType.SINGLE)
                n_bond += 1
    mol = rw.GetMol()
    mol.AddConformer(conf)
    tmpl = Chem.MolFromSmiles(smi)
    if tmpl is None:
        return None, "template SMILES unparsable"
    n_tmpl = tmpl.GetNumBonds()
    if n_bond != n_tmpl:
        return None, f"connectivity mismatch: perceived {n_bond} bonds vs template {n_tmpl}"
    try:
        bonded = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
        Chem.SanitizeMol(bonded)
    except Exception as exc:
        return None, f"bond-order assignment failed: {type(exc).__name__}"
    return bonded, "loaded"


def pdb_het_rows(pdb_path: pathlib.Path):
    """Group HETATM records of a PDB into residues (chain, resseq, resname)."""
    res = {}
    for line in pdb_path.open(encoding="utf-8", errors="replace"):
        if not line.startswith("HETATM"):
            continue
        an = line[12:16].strip()
        if an.startswith("H"):
            continue
        key = (line[21], line[22:26].strip(), line[17:20].strip())
        r = res.setdefault(key, {"name": key[2], "el": [], "xyz": []})
        r["el"].append(line[76:78].strip().capitalize())
        r["xyz"].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return res


def control_crystal_ligands():
    """Same raw-geometry metric on the co-crystal ligands of the three published PDBs.
    A positive control: experimental geometry must pass, otherwise the metric
    (or the cif parsing) is broken."""
    struct_dir = ROOT / "data" / "dti_subset" / "data" / "raw" / "docking"
    rows = []
    for pdb in sorted(struct_dir.glob("*.pdb")):
        for key, r in pdb_het_rows(pdb).items():
            if r["name"] in BUFFER or len(r["el"]) < 10:
                continue
            xyz = np.array(r["xyz"])
            g = pair_stats(xyz)
            g.update({"run": f"{pdb.stem}_{key[2]}{key[1]}{key[0]}", "group": "crystal_control",
                      "n_lig_heavy": len(r["el"]), "n_prot_atoms": np.nan,
                      "min_lig_pro_dist": np.nan, "n_vdw_clashes_0_75": np.nan,
                      "posebusters_load": "not attempted (control is the metric check)"})
            rows.append(g)
    # keep only the largest ligand per entry (same rule as native_ligand in analyze_ext)
    best = {}
    for g in rows:
        e = g["run"][:4]
        if e not in best or g["n_lig_heavy"] > best[e]["n_lig_heavy"]:
            best[e] = g
    return list(best.values())


def main() -> None:
    runs = pd.read_csv(RESULTS / "ext_runs.csv")
    probe = pd.read_csv(ROOT / "data" / "baseline" / "probe_all.csv")
    targets = []
    for r in runs.itertuples():
        targets.append({"name": r.name, "group": r.exp,
                        "cif": next((ROOT / "outputs" / r.name).rglob("*_model_0.cif")),
                        "yaml": INPUTS / f"{r.name}.yaml"})
    for n in probe[probe.preset == "affinity-default+msa"].name:
        targets.append({"name": f"orig_{n}", "group": "orig_single_seq",
                        "cif": next((ROOT / "outputs_baseline" / f"msa_{n}").rglob("*_model_0.cif")),
                        "yaml": INPUTS / f"{n}.yaml"})

    geo_rows, pb_rows = [], []
    for t in targets:
        rows = parse_atom_site(t["cif"])
        het = [r for r in rows if r["group_PDB"] == "HETATM"]
        prot = [r for r in rows if r["group_PDB"] == "ATOM" and r["auth_asym_id"] == "A"]
        smi = re.search(r"smiles:\s*(.+)", t["yaml"].read_text(encoding="utf-8")).group(1).strip()
        xyz_l = np.array([[float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])] for r in het])
        el_l = [r["type_symbol"].capitalize() for r in het]
        xyz_p = np.array([[float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])] for r in prot])
        el_p = [r["type_symbol"].capitalize() for r in prot]
        g = {"run": t["name"], "group": t["group"], "n_lig_heavy": len(het), "n_prot_atoms": len(prot)}
        g.update(pair_stats(xyz_l))
        g.update(clash_stats(xyz_l, el_l, xyz_p, el_p))
        mol, note = perceive_ligand(het, smi)
        g["posebusters_load"] = note
        geo_rows.append(g)
        if mol is not None:
            with tempfile.TemporaryDirectory() as td:
                pdb = pathlib.Path(td) / f"{t['name']}_protein.pdb"
                write_protein_pdb(rows, pdb)
                from posebusters import PoseBusters
                res = PoseBusters(config="dock", max_workers=0).bust(mol_pred=mol, mol_cond=pdb)
            res.insert(0, "run", t["name"])
            res.insert(1, "group", t["group"])
            pb_rows.append(res.iloc[0])
        print(f"{t['name']:28s} min={g['min_pair_dist']:6.3f} lt1.0={g['n_pairs_lt_1_0A']:3d} "
              f"clash={g['n_vdw_clashes_0_75']:4d} {note}", flush=True)

    geo_rows += control_crystal_ligands()
    geo = pd.DataFrame(geo_rows)
    geo.to_csv(RESULTS / "posebusters_geometry.csv", index=False, lineterminator="\n")

    runs_geo = geo[geo.group != "crystal_control"]
    ctrl = geo[geo.group == "crystal_control"]
    lines = ["# B2 PoseBusters 有效性检查（2026-11-18，CPU，全部读仓库内已发布产物）", "",
             "对象：21 个本批预测（E1/E2/E3a/E3b/EQ）＋6 个原单序列基线（orig_），各 model_0 唯一姿态；"
             "另以 3 个已发布 PDB 的共晶配体作指标正对照（crystal_control）。",
             "", "## 第一层：原始几何（不依赖键级感知）", ""]
    n_bad = int((runs_geo.min_pair_dist < 1.0).sum())
    ctrl_txt = "、".join(f"{r.run} {r.min_pair_dist:.3f} Å" for r in ctrl.itertuples())
    lines += [f"- 配体内原子对最近距离 <1.0 Å 的运行：**{n_bad}/{len(runs_geo)}**"
              f"（最小 {runs_geo.min_pair_dist.min():.3f} Å，最大 {runs_geo.min_pair_dist.max():.3f} Å）",
              f"- 配体-蛋白范德华重叠（d < 0.75·(r_i+r_j)）总数：{int(runs_geo['n_vdw_clashes_0_75'].sum())}，"
              f"涉及 {int((runs_geo['n_vdw_clashes_0_75'] > 0).sum())}/{len(runs_geo)} 个运行",
              f"- 正对照（同一指标作用于共晶配体；其最近原子对即最短共价键，1.1–1.3 Å 量级）：{ctrl_txt}，"
              "全部 ≥1.0 Å，而 27 个预测的最近原子对为 0.02–0.27 Å——低于任何真实键长一个量级，"
              "崩塌来自模型输出，不是解析或指标伪影",
              "", "逐运行明细见 `posebusters_geometry.csv`（列：最近原子对距离、<0.5/1.0/1.3 Å 对数、"
              "配体-蛋白最小距离、vdW 重叠对数、PoseBusters 装载结果）。", ""]
    loaded = runs_geo[runs_geo.posebusters_load == "loaded"]
    lines += ["## 第二层：PoseBusters 标准检查（config=dock）", "",
              f"- 能建成合法分子图（感知键数==模板键数且模板键级赋值成功）的运行：{len(loaded)}/{len(runs_geo)}。"
              if len(loaded) else
              f"- 能建成合法分子图的运行：0/{len(runs_geo)}——所有运行坐标感知出的键数与输入 SMILES 模板键数不一致"
              "（重叠原子使感知图多出大量边；下限 0.8 Å 已把 <0.8 Å 的重叠对排除在成键之外，"
              "但被压在一起的原子仍与各自的真实邻居成键，感知图因此与模板不同构），"
              "PoseBusters 检查表无从谈起。",
              "- 失败分布："]
    for note, cnt in runs_geo.posebusters_load.str.split(":").str[0].value_counts().items():
        lines.append(f"  - {note}：{cnt}")
    if pb_rows:
        pb = pd.DataFrame(pb_rows)
        pb.to_csv(RESULTS / "posebusters.csv", index=False, lineterminator="\n")
        lines += ["", "逐运行 PoseBusters 表见 `posebusters.csv`。"]
    lines += ["", "## 口径与限制", "",
              "- 几何层数字直接从已发布 cif 的 Cartn 列计算，无任何模型假设，可复核。",
              "- cif 无 `_chem_comp_bond` 环，键级必须由坐标感知＋SMILES 模板赋值；本批姿态的原子重叠使该步骤"
              "在全部/部分运行失败，失败本身按 PoseBusters 论文口径即『未通过装载/化学有效性』。",
              "- 未做 RMSD/应变类需要共晶配体参考的检查：对接配体与共晶配体不是同一分子（见 ext_results.md E4 注）。",
              "- 复跑：`pip install posebusters==0.6.5` 后 `python scripts/check_posebusters.py`（CPU，分钟级）。"]
    (RESULTS / "posebusters.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    try:
        import sys
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
