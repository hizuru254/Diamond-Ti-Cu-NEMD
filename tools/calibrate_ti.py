#!/usr/bin/env python3
"""Calibrate hcp Ti with NEP89: a and c/a."""
import re, subprocess, sys, json
from pathlib import Path
import numpy as np

BASE = [(0,0,0), (2/3,1/3,0.5)]
M = 47.867
EXE = "/home/cc/lammps/src/lmp_mpi"
MODEL = "/home/cc/si_ti_cu_nemd/potentials/nep89.txt"
W = Path("/tmp/ti_calib"); W.mkdir(parents=True, exist_ok=True)

def energy(a, covera, nx=4, ny=4, nz=3):
    c = a*covera
    # orthorhombic hcp cell: x = a*nx, y = sqrt(3)*a*ny, z = c*nz
    lx, ly, lz = a*nx, a*np.sqrt(3)*ny, c*nz
    pts = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for (bx, by, bz) in BASE:
                    pts.append(((i+bx)*a, (j+by)*a*np.sqrt(3), (k+bz)*c))
    lines = ["Ti hcp", "", f"{len(pts)} atoms", "1 atom types", "",
             f"0 {lx:.6f} xlo xhi", f"0 {ly:.6f} ylo yhi", f"0 {lz:.6f} zlo zhi", "",
             "Masses", "", f"1 {M}", "", "Atoms # atomic", ""]
    for i, p in enumerate(pts, 1):
        lines.append(f"{i} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    (W/"cell.data").write_text("\n".join(lines)+"\n", encoding="utf-8")
    inp = f"""units metal
atom_style atomic
boundary p p p
read_data cell.data
pair_style nep
pair_coeff * * {MODEL} Ti
neighbor 2.0 bin
thermo_style custom step pe
minimize 1e-10 1e-10 3000 30000
"""
    (W/"in.eos").write_text(inp, encoding="utf-8")
    subprocess.run([EXE,"-in","in.eos","-log","log.eos","-screen","none"], cwd=W,
                   check=True, stdout=subprocess.DEVNULL)
    txt = (W/"log.eos").read_text(errors="replace")
    vals = re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", txt, re.M)
    return float(vals[-1])/len(pts)

best = None
for a in np.arange(2.85, 3.05, 0.02):
    for cr in (1.55, 1.59, 1.63):
        e = energy(float(a), cr)
        if best is None or e < best[0]:
            best = (e, float(a), cr)
print(f"coarse: E={best[0]:.4f} a={best[1]:.3f} c/a={best[2]}")
e0, a0, cr0 = best
best = None
for a in np.arange(a0-0.03, a0+0.031, 0.01):
    for cr in np.arange(cr0-0.05, cr0+0.051, 0.02):
        e = energy(float(a), float(cr))
        if best is None or e < best[0]:
            best = (e, float(a), float(cr))
print(f"fine:   E={best[0]:.4f} a={best[1]:.4f} c/a={best[2]:.3f} c={best[1]*best[2]:.4f}")
json.dump({"Ti": {"a0_A": round(best[1],4), "covera": round(best[2],4),
                  "c_A": round(best[1]*best[2],4), "E0_eV_atom": round(best[0],4)}},
          open("/home/cc/si_ti_cu_nemd/potentials/lattice_Ti.json","w"), indent=2)
