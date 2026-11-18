"""Sequential runner for the Boltz-2 extension experiments.

One subprocess per complex (same pattern as the original probe runner): wall
clock and peak GPU memory sampled from nvidia-smi in a side thread, full child
transcript kept, results appended to results/ext_runs.csv after EVERY run so a
crash or manual stop never loses completed rows. Resumable: a run whose
affinity JSON already exists is skipped unless --force.

Presets (structure arm only; affinity head stays at library default 5 x 200):
  E1 / E2 / default : 1 diffusion x 200 steps x 3 recycling (library default)
  E3a               : 1 x 25 x 1
  E3b               : 1 x 100 x 2
Seed 42 everywhere, matching the original probe runs.

Differences vs the local working copy (only path/constants, no protocol):
  - boltz interpreter: env BOLTZ_PY (default: the python you run this with),
    CLI wrapper boltz_cli.py vendored next to this script;
  - model weight cache: env BOLTZ_CACHE (default boltz_cache/ under the repo,
    git-ignored; ~13 GB, downloaded by boltz on first use);
  - inputs/, outputs/, results/ are repo-root directories.
Requires a CUDA GPU (observed peak ~11 GiB at the default preset on an
RTX 2080 Ti; ~14 min per default-preset BRAF complex, 205 s single-sequence).
"""
from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys
import threading
import time

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
VENV_PY = pathlib.Path(os.environ.get("BOLTZ_PY", sys.executable))
BOLTZ_CLI = HERE / "boltz_cli.py"
CACHE = pathlib.Path(os.environ.get("BOLTZ_CACHE", str(ROOT / "boltz_cache")))
RUNS_CSV = ROOT / "results" / "ext_runs.csv"

PRESETS = {
    "E3a": ["--diffusion_samples", "1", "--recycling_steps", "1", "--sampling_steps", "25"],
    "E3b": ["--diffusion_samples", "1", "--recycling_steps", "2", "--sampling_steps", "100"],
}


class Sampler(threading.Thread):
    def __init__(self) -> None:
        super().__init__(daemon=True)
        self.stop_flag = threading.Event()
        self.peak_mib = 0

    def run(self) -> None:
        proc = subprocess.Popen(
            ["nvidia-smi", "--query-gpu=memory.used,utilization.gpu",
             "--format=csv,noheader,nounits", "-lms", "500"],
            stdout=subprocess.PIPE, text=True, bufsize=1,
        )
        for line in proc.stdout:
            if self.stop_flag.is_set():
                break
            try:
                used, _util = (float(x) for x in line.split(","))
            except ValueError:
                continue
            self.peak_mib = max(self.peak_mib, used)
        proc.terminate()

    def stop(self) -> None:
        self.stop_flag.set()


def run_one(yaml_path: pathlib.Path, out_dir: pathlib.Path, extra: list[str]):
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(VENV_PY), str(BOLTZ_CLI), "predict", str(yaml_path),
           "--cache", str(CACHE), "--out_dir", str(out_dir),
           "--accelerator", "gpu", "--devices", "1", "--seed", "42", "--override", *extra]
    s = Sampler()
    s.start()
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    s.stop_flag.set()
    s.join(timeout=5)
    aff, conf = {}, {}
    for f in out_dir.rglob("affinity_*.json"):
        aff = json.loads(f.read_text(encoding="utf-8"))
    for f in out_dir.rglob("confidence_*.json"):
        conf = json.loads(f.read_text(encoding="utf-8"))
    (out_dir / "child.log").write_text(
        f"# exit={proc.returncode}\n# cmd={' '.join(cmd)}\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}",
        encoding="utf-8")
    if proc.returncode != 0:
        print(f"!! exit {proc.returncode}\n{(proc.stdout + proc.stderr)[-1200:]}", flush=True)
    return dt, s.peak_mib, aff, conf


def main() -> None:
    force = "--force" in sys.argv
    only = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--only=")), None)
    meta = pd.read_csv(ROOT / "results" / "ext_meta.csv").set_index("name")
    done = set()
    if RUNS_CSV.exists() and not force:
        done = set(pd.read_csv(RUNS_CSV)["name"])
    names = sorted(meta.index)
    if only:
        names = [n for n in names if only in n]
    for name in names:
        m = meta.loc[name]
        out_dir = ROOT / "outputs" / name
        if name in done and list(out_dir.rglob("affinity_*.json")):
            print(f"== {name}: already done, skip", flush=True)
            continue
        extra = PRESETS.get(m["exp"], [])
        if m["msa_mode"] == "server":
            extra = [*extra, "--use_msa_server"]
        print(f">> {name} ({m['exp']}, msa={m['msa_mode']})", flush=True)
        dt, peak, aff, conf = run_one(ROOT / "inputs" / f"{name}.yaml", out_dir, extra)
        row = {
            "name": name, "exp": m["exp"], "gene": m["gene"], "pdb_id": m["pdb_id"],
            "construct_len": m["construct_len"], "ligand_id": m["ligand_id"],
            "true_pkd": m["true_pkd"], "vina_best": m["vina_best"], "msa_mode": m["msa_mode"],
            "structure_preset": m["exp"] if m["exp"] in PRESETS else "default",
            "seconds": round(dt, 1), "peak_gpu_MiB": int(peak),
            "affinity_pred": aff.get("affinity_pred_value"),
            "binary_bind_prob": aff.get("affinity_probability_binary"),
            "ptm": conf.get("ptm"), "iptm": conf.get("iptm"),
            "complex_plddt": conf.get("complex_plddt"), "complex_ipde": conf.get("complex_ipde"),
            "confidence_score": conf.get("confidence_score"),
            "affinity_ens1": aff.get("affinity_pred_value1"),
            "affinity_ens2": aff.get("affinity_pred_value2"),
        }
        RUNS_CSV.parent.mkdir(exist_ok=True)
        pd.DataFrame([row]).to_csv(RUNS_CSV, mode="a", header=not RUNS_CSV.exists(), index=False)
        print(f"   {dt:.0f}s peak {peak:.0f} MiB affinity={aff.get('affinity_pred_value')} "
              f"ptm={conf.get('ptm')} iptm={conf.get('iptm')}", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
