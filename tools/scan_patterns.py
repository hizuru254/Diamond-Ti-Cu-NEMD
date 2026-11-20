#!/usr/bin/env python3
"""Unattended pattern scan: build -> gate -> run -> measure -> rank, sequentially.

Runs are strictly serial because the NEP pair style is memory-bandwidth bound:
two concurrent jobs drop to 6-14 steps/s each instead of 61 for a single job.

Usage:
  python3 tools/scan_patterns.py --designs "0,0,0:FLAT" "6,6,10.325:IL_H6W6P10" ...
"""
from __future__ import annotations
import argparse, json, re, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable


def sh(cmd: str, cwd: Path = ROOT, timeout: int = 7200) -> tuple[int, str]:
    r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    return r.returncode, (r.stdout or "") + (r.stderr or "")


def build_and_gate(h: float, w: float, p: float, tag: str) -> bool:
    cmd = (f"{PY} tools/build_pair_hex.py --interlock {h} --groove-w {w} --groove-p {p} "
           f"--out-structure structures/{tag}.data --out-geom structures/{tag}.geom.json")
    rc, out = sh(cmd)
    print(out.strip().splitlines()[-1] if out.strip() else "", flush=True)
    if "PASS" not in out:
        print(f"[SCAN] {tag}: geometry check failed"); return False
    rc, out = sh(f"{PY} tools/preflight.py --config config/params.json "
                 f"--geom structures/{tag}.geom.json --structure structures/{tag}.data "
                 f"--potentials potentials/potentials.in --dry-run")
    ok = "READY" in out
    print(f"[SCAN] {tag}: preflight {'READY' if ok else 'BLOCKED'}", flush=True)
    if not ok:
        print(out[-800:])
    return ok


def make_run(tag: str) -> str | None:
    rc, out = sh(f"grep -c '^{tag},' scan.csv || true")
    if rc != 0 or out.strip() in ("0", ""):
        with open(ROOT/"scan.csv", "a", encoding="utf-8") as f:
            f.write(f"{tag},pair,0,0,0,0,0,single,pending,scan design {tag}\n")
    rc, out = sh(f"{PY} tools/new_run.py --config config/params.json --scan-id {tag} "
                 f"--geom structures/{tag}.geom.json --structure structures/{tag}.data "
                 f"--potentials potentials/potentials.in --scan-csv scan.csv")
    m = re.search(r"RUN_DIR=(\S+)", out)
    return m.group(1) if m else None


def run_one(run_dir: str) -> bool:
    rc, out = sh(f'MPIRUN_EXTRA="--cpu-set 0-23" bash tools/run_seed.sh {run_dir} 1001')
    return "exit=0" in out


def measure(run_dir: str) -> dict:
    rc, out = sh(f"{PY} tools/tbr_extrap.py --run-dir {run_dir} --seed seed_1001", timeout=1800)
    m = re.search(r"interface TBR \(z->0\) = ([0-9.eE+-]+)", out)
    r2 = re.search(r"R2=([0-9.]+)", out)
    return {"tbr": float(m.group(1)) if m else None,
            "r2": float(r2.group(1)) if r2 else None}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", nargs="+", required=True, help="H,W,P:tag")
    ap.add_argument("--out", default="reports/pattern_scan.json")
    a = ap.parse_args()
    results = {}
    for d in a.designs:
        spec, tag = d.split(":")
        h, w, p = [float(x) for x in spec.split(",")]
        print(f"\n===== {tag}  H={h} W={w} P={p} =====", flush=True)
        if not build_and_gate(h, w, p, tag):
            results[tag] = {"error": "gate failed"}; continue
        run_dir = make_run(tag)
        if not run_dir:
            results[tag] = {"error": "new_run failed"}; continue
        t0 = time.time()
        ok = run_one(run_dir)
        dt = time.time() - t0
        res = measure(run_dir) if ok else {"error": "run failed"}
        res["wall_s"] = round(dt)
        res["run_dir"] = run_dir
        results[tag] = res
        print(f"[SCAN] {tag}: {res}", flush=True)
        Path(ROOT/a.out).write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("\n=== ranking ===")
    ok = {k: v for k, v in results.items() if v.get("tbr")}
    for k, v in sorted(ok.items(), key=lambda kv: kv[1]["tbr"]):
        print(f"  {k:16s} TBR={v['tbr']:.3e}  R2={v['r2']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
