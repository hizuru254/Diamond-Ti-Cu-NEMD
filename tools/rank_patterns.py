#!/usr/bin/env python3
"""Collect every measured run into one ranked table (TBR, R2, stability, atoms)."""
from __future__ import annotations
import json, re, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable


def measure(run: Path) -> dict:
    out = {"run": run.name}
    r = subprocess.run([PY, "tools/tbr_extrap.py", "--run-dir", str(run),
                        "--seed", "seed_1001"], cwd=ROOT, capture_output=True, text=True)
    m = re.search(r"interface TBR \(z->0\) = ([0-9.eE+-]+)", r.stdout)
    r2 = re.search(r"R2=([0-9.]+)", r.stdout)
    out["tbr"] = float(m.group(1)) if m else None
    out["r2"] = float(r2.group(1)) if r2 else None
    sd = run/"seed_1001"
    log = sd/"log.pattern"
    if log.exists():
        pe = []
        for line in log.read_text(errors="replace").splitlines():
            s = line.split()
            if len(s) == 7 and s[0].isdigit():
                try: pe.append((int(s[0]), float(s[2])))
                except ValueError: pass
        if pe:
            a = np.array(pe); half = a[len(a)//2:]
            out["dpe_dt"] = float(np.polyfit(half[:,0], half[:,1], 1)[0]*1000)
    res = sd/"result.json"
    if res.exists():
        try:
            d = json.loads(res.read_text())
            out["drift_K"] = d.get("temperature", {}).get("drift_K")
        except Exception:
            pass
    meta = run/"run_meta.json"
    if meta.exists():
        try:
            d = json.loads(meta.read_text())
            out["atoms"] = d.get("inputs", {}).get("structure_atoms")
        except Exception:
            pass
    return out


def main() -> int:
    rows = []
    for run in sorted((ROOT/"workspace").glob("run_v*")):
        if not (run/"seed_1001"/"log.pattern").exists():
            continue
        rows.append(measure(run))
    rows = [r for r in rows if r.get("tbr")]
    rows.sort(key=lambda r: r["tbr"])
    base = min(r["tbr"] for r in rows)
    print(f"{'run':10s} {'TBR (m2K/W)':>13s} {'rel':>6s} {'R2':>6s} {'dPE/dt':>9s}")
    for r in rows:
        print(f"{r['run']:10s} {r['tbr']:13.3e} {r['tbr']/base:6.2f} "
              f"{r.get('r2') or 0:6.3f} {r.get('dpe_dt') or 0:+9.4f}")
    out = ROOT/"reports/pattern_ranking.json"
    out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"\n[OK] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
