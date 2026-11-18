#!/usr/bin/env python3
"""Bulk calibration for the Cu/X/diamond NEP (elements: Cu C Cr Mo Ti W).

Prints the FULL energy surface so a non-convex or edge-hugging minimum is visible
instead of being silently accepted.
"""
from __future__ import annotations
import re, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path("/home/cc/si_ti_cu_nemd")
EXE = "/home/cc/lammps/src/lmp_mpi"
MODEL = ROOT / "potentials/CuTiC_NEP.txt"
W = Path("/tmp/calib2"); W.mkdir(parents=True, exist_ok=True)

BASIS = {
    "C":  [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0),(.25,.25,.25),(.25,.75,.75),(.75,.25,.75),(.75,.75,.25)],
    "Cu": [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)],
}
MASS = {"C": 12.011, "Cu": 63.546, "Ti": 47.867}


def write_cell(el, a, c=None, covera=None, n=3):
    if el == "Ti":
        c = c if c is not None else a * covera
        basis = [(0,0,0), (2/3,1/3,.5)]
        lx, ly, lz = a*n, a*np.sqrt(3)*n, c*n
        pts = [((i+bx)*a, (j+by)*a*np.sqrt(3), (k+bz)*c)
               for i in range(n) for j in range(n) for k in range(n) for bx,by,bz in basis]
    else:
        basis = BASIS[el]
        lx = ly = lz = a*n
        pts = [((i+bx)*a, (j+by)*a, (k+bz)*a)
               for i in range(n) for j in range(n) for k in range(n) for bx,by,bz in basis]
    lines = [el, "", f"{len(pts)} atoms", "1 atom types", "",
             f"0 {lx:.6f} xlo xhi", f"0 {ly:.6f} ylo yhi", f"0 {lz:.6f} zlo zhi", "",
             "Masses", "", f"1 {MASS[el]}", "", "Atoms # atomic", ""]
    for i, p in enumerate(pts, 1):
        lines.append(f"{i} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    (W/"cell.data").write_text("\n".join(lines)+"\n", encoding="utf-8")
    return len(pts)


def energy(el, a, c=None, covera=None):
    n = write_cell(el, a, c, covera)
    inp = f"""units metal
atom_style atomic
boundary p p p
read_data cell.data
pair_style nep
pair_coeff * * {MODEL} {el}
neighbor 2.0 bin
neigh_modify one 4000
thermo_style custom step pe
minimize 1e-12 1e-12 5000 50000
print "EPER $(pe)"
"""
    (W/"in.e").write_text(inp, encoding="utf-8")
    r = subprocess.run([EXE, "-in", "in.e", "-log", "log.e", "-screen", "none"],
                       cwd=W, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    txt = (W/"log.e").read_text(errors="replace")
    m = re.findall(r"EPER\s+(-?\d+\.\d+)", txt)
    if not m:
        return None
    return float(m[-1]) / n


def scan_cubic(el, lo, hi, step):
    print(f"--- {el} ---")
    best = None
    for a in np.arange(lo, hi + 1e-9, step):
        e = energy(el, float(a))
        if e is None:
            print(f"  a={a:.3f}  FAILED"); continue
        print(f"  a={a:.3f}  E={e:.4f}")
        if best is None or e < best[1]:
            best = (float(a), e)
    return best


print("=== Cu (fcc), exp a=3.615, E0=-3.49 ===")
cu = scan_cubic("Cu", 3.40, 3.90, 0.05)
print("=== C (diamond), exp a=3.567, E0=-7.37 ===")
c = scan_cubic("C", 3.40, 3.80, 0.05)
print("=== Ti (hcp) 2-D: exp a=2.951, c/a=1.587, E0=-4.85 ===")
best = None
for a in np.arange(2.85, 3.10, 0.05):
    row = []
    for cr in (1.50, 1.55, 1.587, 1.63):
        e = energy("Ti", float(a), covera=float(cr))
        row.append(e)
        if e is not None and (best is None or e < best[0]):
            best = (e, float(a), float(cr))
    print(f"  a={a:.2f}  " + "  ".join(f"c/a={cr}:{('%.3f'%v) if v is not None else 'FAIL'}"
                                       for cr, v in zip((1.50,1.55,1.587,1.63), row)))
print(f"[MIN] Cu a0={cu[0]:.3f} E={cu[1]:.4f} | C a0={c[0]:.3f} E={c[1]:.4f} | "
      f"Ti a0={best[1]:.3f} c/a={best[2]:.3f} E={best[0]:.4f}")
