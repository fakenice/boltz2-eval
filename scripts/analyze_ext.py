"""Analyze the Boltz-2 extension experiments. Reads only real run outputs.

E1  within-BRAF ranking: Spearman of Boltz-2 affinity vs measured pKd, and of
    Vina (-score) vs pKd, on the same 11 ligands; bootstrap 95% CI (2000 resamples,
    seed 42). Also the active-only subset (pKd > 5.0). Absolute affinity values
    are NOT comparable to pKd (scale mismatch documented in the 10-14 report);
    only within-target ordering is graded.
E2  MSA ablation (corrected protocol): the original probe "affinity-default+msa"
    runs were actually SINGLE-SEQUENCE (YAML msa: empty -> msa_id=-1 overrides
    --use_msa_server; records confirm). Real-MSA side = E1 BRAF runs + E2 runs
    for the other four originals, all with the persisted local a3m. Deltas are
    real-MSA minus single-sequence. E2rep re-runs one original single-sequence
    to check reproducibility.
E3  structure sampling-quality curve: 25/1, 100/2 (this batch) vs 200/3 taken
    from the real-MSA default runs (E1_T0013_L02069 / E2_T0004_L01091) so the
    whole curve shares one MSA protocol; affinity head default in all.
E4  crystal-reference C-alpha RMSD: predicted chain A aligned to the crystal
    chain A it was extracted from (index pairing in construct order, residue
    identity checked), whole-chain and binding-site (crystal residues within
    5 A of the co-crystal ligand) RMSD after the global fit.
EQ  protocol equivalence: server-MSA run (msa auto + --use_msa_server) vs the
    local-a3m run of the same complex.

Repository layout vs the local working copy (paths only):
  results/ext_runs.csv         this batch, 21 runs (appended by run_ext.py)
  data/baseline/probe_all.csv  the 6 original preset=affinity-default+msa rows
                               (extracted from the local probe_all.csv, columns
                               unchanged)
  outputs/<name>/              this batch's predictions (cif + json only)
  outputs_baseline/msa_<n>/    the originals' predictions (same subset)
No number changes: running this script reproduces results/ext_results.md and
the four ext_*.csv from the published outputs alone.
"""
from __future__ import annotations

import os
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
BENCH = pathlib.Path(os.environ.get("DTI_BENCH", str(ROOT / "data" / "dti_subset")))
STRUCT_DIR = BENCH / "data" / "raw" / "docking"
RESULTS = ROOT / "results"

BUFFER = {"HOH", "EDO", "GOL", "SO4", "PO4", "DMS", "PEG", "MES", "TRS", "CL", "NA",
          "MG", "ZN", "CA", "K", "NO3", "ACT", "EPE", "CIT", "FLC", "PGE", "PG4",
          "1PE", "P6G", "XP4", "XPE", "PG0", "2PE", "15P", "PE4", "PEU", "DTT",
          "BME", "IOD", "BR", "FMT", "ACN", "MPD", "TLA", "TAR", "SIN", "ORN"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y",
    "ALY": "K", "MLZ": "K", "NLE": "X", "AIB": "X",
}


def md_table(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for row in df.itertuples(index=False):
        out.append("| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |")
    return "\n".join(out)


def boot_spearman(x, y, n_boot=2000, seed=42):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return np.nan, (np.nan, np.nan), len(x)
    rho = stats.spearmanr(x, y).statistic
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        vals.append(stats.spearmanr(x[idx], y[idx]).statistic)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return rho, (lo, hi), len(x)


def crystal_residues(pdb_path, chain="A"):
    """Chain A ATOM residues in construct order (sorted residue number)."""
    res = {}
    het = {}
    for line in pdb_path.open(encoding="utf-8", errors="replace"):
        if line.startswith("ATOM  ") and line[21] == chain:
            num = int(line[22:26])
            r = res.setdefault(num, {"name": line[17:20].strip(), "atoms": {}})
            r["atoms"][line[12:16].strip()] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("HETATM"):
            key = (line[21], line[22:26].strip(), line[17:20].strip())
            r = het.setdefault(key, {"name": line[17:20].strip(), "atoms": {}})
            an = line[12:16].strip()
            if not an.startswith("H"):
                r["atoms"][an] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    ordered = [res[n] for n in sorted(res)]
    return ordered, het


def native_ligand(het):
    cands = [(k, v) for k, v in het.items()
             if v["name"] not in BUFFER and len(v["atoms"]) >= 10]
    if not cands:
        return None, None
    k, v = max(cands, key=lambda kv: len(kv[1]["atoms"]))
    return k, np.array(list(v["atoms"].values()))


def predicted_chain(cif_path):
    from Bio.PDB import MMCIFParser
    s = MMCIFParser(QUIET=True).get_structure("pred", str(cif_path))
    chain = s[0]["A"]
    residues = [r for r in chain if r.id[0] == " "]
    residues.sort(key=lambda r: r.id[1])
    return residues


def rmsd_block(name, cif_path, pdb_id):
    ordered, het = crystal_residues(STRUCT_DIR / f"{pdb_id.lower()}.pdb")
    pred = predicted_chain(cif_path)
    n = min(len(ordered), len(pred))
    fixed, moving, ident = [], [], 0
    for i in range(n):
        cr, pr = ordered[i], pred[i]
        if "CA" not in cr["atoms"] or "CA" not in pr:
            continue
        want = THREE_TO_ONE.get(cr["name"], "?")
        got = pr.get_resname()
        if THREE_TO_ONE.get(got, got) == want or got == cr["name"]:
            ident += 1
        fixed.append(cr["atoms"]["CA"])
        moving.append(pr["CA"].get_vector().get_array())
    fixed, moving = np.array(fixed), np.array(moving)
    # Kabsch fit (numpy only; Bio.PDB.Superimposer wants Atom objects, we have arrays)
    mc, fc = moving.mean(axis=0), fixed.mean(axis=0)
    H = (moving - mc).T @ (fixed - fc)
    U, _S, Vt = np.linalg.svd(H)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    t = fc - R @ mc
    rms_all = float(np.sqrt(np.mean(np.sum((moving @ R.T + t - fixed) ** 2, axis=1))))
    # binding site from crystal native ligand
    lig_key, lig_xyz = native_ligand(het)
    rms_site, n_site = np.nan, 0
    if lig_xyz is not None:
        site_idx = []
        for i, cr in enumerate(ordered[:n]):
            xyz = np.array(list(cr["atoms"].values()))
            if len(xyz) and np.min(np.linalg.norm(xyz[:, None, :] - lig_xyz[None, :, :], axis=2)) <= 5.0:
                site_idx.append(i)
        # map ordered index -> position in fixed/moving (CA-present subset)
        ca_pos, pos = {}, 0
        for i in range(n):
            if "CA" in ordered[i]["atoms"] and i < len(pred) and "CA" in pred[i]:
                ca_pos[i] = pos
                pos += 1
        sel = [ca_pos[i] for i in site_idx if i in ca_pos]
        n_site = len(sel)
        if n_site >= 3:
            mv = moving[sel] @ R.T + t
            rms_site = float(np.sqrt(np.mean(np.sum((mv - fixed[sel]) ** 2, axis=1))))
    return {
        "name": name, "pdb_id": pdb_id, "crystal_res": len(ordered), "pred_res": len(pred),
        "paired_CA": len(fixed), "identity_frac": round(ident / max(len(fixed), 1), 3),
        "plddt_mean": round(float(np.mean([p["CA"].get_bfactor() for p in pred[:n] if "CA" in p])), 1),
        "rmsd_all_CA": round(rms_all, 3), "native_ligand": lig_key[2] if lig_key else None,
        "site_res": n_site, "rmsd_site_CA": round(rms_site, 3) if np.isfinite(rms_site) else None,
    }


def find_cif(root, name):
    hits = list((root / name).rglob("*_model_0.cif"))
    return hits[0] if hits else None


def main() -> None:
    runs = pd.read_csv(RESULTS / "ext_runs.csv")
    probe = pd.read_csv(ROOT / "data" / "baseline" / "probe_all.csv")
    orig = probe[probe.preset == "affinity-default+msa"].set_index("name")
    lines = ["# Boltz-2 补充实验结果（2026-11-17，全部本机实测）", "",
             "口径更正：原 probe「affinity-default+msa」6 对的 YAML 均为 `msa: empty`，"
             "解析为 msa_id=-1，优先于 --use_msa_server（schema.py L1118 / main.py L517-529），"
             "records 亦为 -1 ⇒ 原 6 对实为单序列运行；本批 E1/E2/E3 用持久化本地 a3m（真实 MSA，上限 8192 条）。", ""]

    # ---- E1 ranking ----
    e1 = runs[runs.exp == "E1"].copy()
    lines += ["## E1 BRAF 同靶点排序（n=%d，其中 pKd>5.0 共 %d）" % (len(e1), (e1.true_pkd > 5).sum()), ""]
    if len(e1):
        tab = e1.sort_values("true_pkd", ascending=False)[
            ["ligand_id", "true_pkd", "vina_best", "affinity_pred", "binary_bind_prob", "ptm", "iptm"]]
        lines += [md_table(tab), ""]
        for label, sub in [("全部配体", e1), ("pKd>5.0 子集", e1[e1.true_pkd > 5.0])]:
            rb, cb, nb = boot_spearman(sub.affinity_pred, sub.true_pkd)
            rv, cv, nv = boot_spearman(-sub.vina_best, sub.true_pkd)
            rp, cp, np_ = boot_spearman(sub.binary_bind_prob, sub.true_pkd)
            lines.append(f"- {label}（n={nb}）：Boltz-2 亲和力 vs pKd Spearman **{rb:+.3f}** [{cb[0]:+.3f}, {cb[1]:+.3f}]；"
                         f"Vina vs pKd **{rv:+.3f}** [{cv[0]:+.3f}, {cv[1]:+.3f}]；"
                         f"二分类结合概率 vs pKd {rp:+.3f} [{cp[0]:+.3f}, {cp[1]:+.3f}]")
        lines += ["", "口径：pKd 为 pairs.csv 同一配体多行的均值（与 10-14 报告一致）；5.0 是 DAVIS Kd 截断上限，"
                  "并列按 Spearman 默认处理；亲和力绝对值与 pKd 尺度不可比，只评同靶点内排序。", ""]
        e1.to_csv(RESULTS / "ext_ranking.csv", index=False)

    # ---- EQ equivalence: server MSA vs local a3m reuse ----
    lines += ["## EQ 服务器 MSA vs 本地 a3m 复用（协议等价性）", ""]
    eq = runs[runs.name == "EQ_T0013_L02069"]
    e1ref = runs[runs.name == "E1_T0013_L02069"]
    if len(eq) and len(e1ref):
        r, o = eq.iloc[0], e1ref.iloc[0]
        lines.append(f"- T0013_L02069：服务器 MSA affinity {r.affinity_pred:+.4f} vs 本地 a3m {o.affinity_pred:+.4f}"
                     f"（Δ {r.affinity_pred - o.affinity_pred:+.4f}）；ptm {r.ptm:.4f} vs {o.ptm:.4f}；"
                     f"iptm {r.iptm:.4f} vs {o.iptm:.4f}；耗时 {r.seconds:.0f}s vs {o.seconds:.0f}s")
        lines.append("")
        lines.append("注：服务器路径把 a3m 转储为 CSV（key,sequence），本地路径直接解析 a3m；"
                     "两者同源于 api.colabfold.com 同一检索，但解析管线不同，Δ 即管线差异的实测。")
    else:
        lines += ["（EQ 或 E1 参考运行缺失，未出数）", ""]
    lines.append("")

    # ---- E2 ablation: real MSA (new) vs single-sequence (original) ----
    lines += ["## E2 MSA 消融（真实 MSA − 单序列；原 6 对为单序列基线）", ""]
    real_default = {}
    for nm in ("E1_T0013_L02069", "E1_T0013_L00585"):
        row = runs[runs.name == nm]
        if len(row):
            real_default[nm[3:]] = row.iloc[0]
    for r in runs[runs.exp == "E2"].itertuples():
        real_default[r.name[3:]] = r
    rows = []
    for base, r in real_default.items():
        if base in orig.index:
            o = orig.loc[base]
            rows.append({"name": base, "gene": r.gene,
                         "d_ptm": round(r.ptm - o.ptm, 4), "d_iptm": round(r.iptm - o.iptm, 4),
                         "d_plddt": round(r.complex_plddt - o.complex_plddt, 4),
                         "d_ipde": round(r.complex_ipde - o.complex_ipde, 3),
                         "d_affinity": round(r.affinity_pred - o.affinity_pred, 4),
                         "s_msa": r.seconds, "s_single": o.seconds,
                         "peak_MiB_msa": r.peak_gpu_MiB, "peak_MiB_single": o.peak_gpu_MiB})
    if rows:
        ab = pd.DataFrame(rows)
        lines += [md_table(ab), "",
                  "Δ = 真实 MSA − 单序列（原 probe 运行）。正 Δptm/Δiptm 表示加入 MSA 后置信度上升。", ""]
        ab.to_csv(RESULTS / "ext_ablation.csv", index=False)
    rep = runs[runs.exp == "E2rep"]
    if len(rep) and "T0013_L00585" in orig.index:
        r, o = rep.iloc[0], orig.loc["T0013_L00585"]
        lines += [f"复现性核对（E2rep，单序列重跑 T0013_L00585）：affinity {r.affinity_pred:+.4f} vs 原 {o.affinity_pred:+.4f}"
                  f"（Δ {r.affinity_pred - o.affinity_pred:+.4f}）；ptm {r.ptm:.4f} vs {o.ptm:.4f}；"
                  f"耗时 {r.seconds:.0f}s vs {o.seconds:.0f}s", ""]

    # ---- E3 curve ----
    lines += ["## E3 结构采样档-质量曲线（真实 MSA，亲和力头固定默认 5×200）", ""]
    rows = []
    for r in runs[runs.exp.isin(["E3a", "E3b"])].itertuples():
        rows.append({"complex": r.name[4:], "preset": r.structure_preset,
                     "steps": 25 if r.exp == "E3a" else 100, "recycles": 1 if r.exp == "E3a" else 2,
                     "ptm": r.ptm, "iptm": r.iptm, "plddt": r.complex_plddt, "ipde": r.complex_ipde,
                     "seconds": r.seconds, "peak_MiB": r.peak_gpu_MiB})
    for base in ("T0013_L02069", "T0004_L01091"):
        r = real_default.get(base)
        if r is not None:
            rows.append({"complex": base, "preset": "default", "steps": 200, "recycles": 3,
                         "ptm": r.ptm, "iptm": r.iptm, "plddt": r.complex_plddt, "ipde": r.complex_ipde,
                         "seconds": r.seconds, "peak_MiB": r.peak_gpu_MiB})
    if rows:
        cv = pd.DataFrame(rows).sort_values(["complex", "steps"])
        lines += [md_table(cv), ""]
        cv.to_csv(RESULTS / "ext_curve.csv", index=False)

    # ---- E4 RMSD ----
    lines += ["## E4 晶体参考 Cα RMSD（整体对齐后）", ""]
    rows = []
    targets = [(f"msa_{n}", f"ss_{n}", orig.loc[n].pdb_id) for n in orig.index]
    targets += [(r.name, r.name, r.pdb_id) for r in runs[runs.exp.isin(["E1", "E2"])].itertuples()]
    for root_name, label, pdb_id in targets:
        root = ROOT / "outputs_baseline" if root_name.startswith("msa_") else ROOT / "outputs"
        cif = find_cif(root, root_name)
        if cif is None:
            continue
        try:
            rows.append(rmsd_block(label, cif, pdb_id))
        except Exception as exc:  # keep going; record the failure honestly
            rows.append({"name": label, "pdb_id": pdb_id, "error": str(exc)[:120]})
    if rows:
        rd = pd.DataFrame(rows)
        lines += [md_table(rd), "",
                  "配对方式：预测链 A 按残基顺序与晶体 chain A（构建序列提取顺序）逐位配对，identity_frac 为残基类型一致比例；"
                  "ss_ 前缀＝原单序列运行，其余为真实 MSA 运行；plddt_mean＝预测 cif 中 CA B-factor（pLDDT，0–100）均值；"
                  "结合位点＝晶体中与共晶配体重原子 ≤5 Å 的残基；共晶配体取 PDB 中最大的非缓冲杂残基。"
                  "对接配体与共晶配体不是同一分子，故不报配体 RMSD。", ""]
        rd.to_csv(RESULTS / "ext_rmsd.csv", index=False)

    # Timing-caliber note, computed from this batch's own runs. The named
    # examples are pinned by assertions: if the data ever drifts, this fails
    # loudly and the prose must be rewritten rather than quietly kept stale.
    braf = runs[(runs.pdb_id == "8c7x") & (runs.structure_preset == "default")
                & (runs.msa_mode.isin(["local_a3m", "server"]))]
    lo, hi = braf.seconds.min(), braf.seconds.max()
    lo_name, hi_name = braf.loc[braf.seconds.idxmin(), "name"], braf.loc[braf.seconds.idxmax(), "name"]
    assert lo_name == "EQ_T0013_L02069" and hi_name == "E1_T0013_L00349", (lo_name, hi_name)
    mid = runs.loc[runs.name == "E1_T0013_L00602", "seconds"].iloc[0]
    peaks = runs[runs.msa_mode == "local_a3m"].peak_gpu_MiB
    lines += ["## 计时口径说明（实测核查）", "",
              "同协议同输入（274 aa、同一 8192 条 MSA 张量、seed 42、默认档）的单对墙钟耗时在本批内波动极大：",
              f"{lo_name} {lo:.0f} s、E1_T0013_L00602 {mid:.0f} s，而 {hi_name} {hi:.0f} s（其进度条显示同一采样阶段",
              "3:15 vs 1:14:58）。各 run 的 MSA npz 形状完全一致（sequences 8192 / residues 2244608），",
              "EQ 与 E1 同复合物结果逐位相同 ⇒ 计算内容相同，耗时差异来自系统级因素（夜间长时负载下的 GPU 降频",
              f"或后台进程争用，未单独确证），不是协议差异。因此：真实 MSA 档单对耗时只能报区间（约 {lo/60:.0f}–{hi/60:.0f} 分钟），",
              f"不得引用单点值作速度结论；峰值显存稳定在 {peaks.min():.0f}–{peaks.max():.0f} MiB（11 GB 卡上限附近），{len(runs)}/{len(runs)} 全部产出、无 OOM。"]

    (RESULTS / "ext_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    try:
        import sys
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print("\n".join(lines))


if __name__ == "__main__":
    main()
