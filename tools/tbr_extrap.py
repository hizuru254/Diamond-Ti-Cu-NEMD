#!/usr/bin/env python3
"""Interface resistance by plane extrapolation.

TBR_eff(z) = dT(+-z)/q includes the bulk resistance 2z/kappa between the two
planes, so plotting TBR_eff against z and extrapolating to z->0 removes the bulk
term and leaves the interface resistance (standard NEMD length extrapolation).
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np


def plane_tbr(run: Path, seed: str, z: float, tool: Path) -> tuple[float, float]:
    r = subprocess.run([sys.executable, str(tool), "--run-dir", str(run),
                        "--seed", seed, "--z", str(z)],
                       capture_output=True, text=True)
    out = r.stdout + r.stderr
    for line in out.splitlines():
        if "TBR_eff" in line:
            try:
                val = float(line.split("TBR_eff =")[1].split("+/-")[0])
                err = float(line.split("+/-")[1].split("(")[0])
                return val, err
            except Exception:
                pass
    return float("nan"), float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--seed", default="seed_1001")
    ap.add_argument("--zs", default="8,12,16,20,24,28,32,36,40")
    ap.add_argument("--tool", default="tools/tbr_plane.py")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    run = Path(a.run_dir)
    zs = [float(x) for x in a.zs.split(",")]
    vals, errs = [], []
    for z in zs:
        v, e = plane_tbr(run, a.seed, z, Path(a.tool))
        vals.append(v); errs.append(e)
        print(f"  z={z:5.1f}  TBR_eff={v:.4e} +/- {e:.1e}")
    vals = np.array(vals); zs = np.array(zs)
    ok = np.isfinite(vals)
    if ok.sum() < 3:
        print("[FAIL] not enough valid points"); return 1
    co = np.polyfit(zs[ok], vals[ok], 1)
    tbr0 = co[1]
    resid = vals[ok] - np.polyval(co, zs[ok])
    r2 = 1 - (resid**2).sum()/((vals[ok]-vals[ok].mean())**2).sum()
    print(f"\n[EXTRAP] TBR(z) = {co[0]:.3e}*z + {co[1]:.3e}   R2={r2:.4f}")
    print(f"[EXTRAP] interface TBR (z->0) = {tbr0:.3e} m2K/W")
    print(f"[EXTRAP] implied bulk kappa = {1/(2*co[0]):.2f} W/mK (both sides averaged)")
    if a.out:
        Path(a.out).write_text(json.dumps({"zs": zs[ok].tolist(), "tbr_eff": vals[ok].tolist(),
                                           "slope": co[0], "intercept": tbr0, "r2": r2},
                                          indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
