#!/usr/bin/env python3
"""hcp Ti calibration with the Cu/Ti/C/diamond NEP: full 35-point grid, full minimize.

Exact triclinic hexagonal cell (a1=(a,0,0), a2=(-a/2,sqrt3 a/2,0), c) with atoms at
(0,0,0) and (1/3,2/3,1/2).  Every point is minimised to 1e-12 and the whole
surface is printed so a non-convex minimum cannot hide.
"""
import re, subprocess, json, time
from pathlib import Path
import numpy as np

W = Path("/tmp/ti_grid"); W.mkdir(parents=True, exist_ok=True)
EXE = "/home/cc/lammps/src/lmp_mpi"
MODEL = "/home/cc/si_ti_cu_nemd/potentials/CuTiC_NEP.txt"
NX, NY, NZ = 3, 3, 2

def energy(a, cr):
    c = a*cr
    lx, ly, lz = NX*a, NY*a*np.sqrt(3)/2, NZ*c
    xy = -a/2
    atoms = []
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                for fx, fy, fz in ((0,0,0), (1/3,2/3,.5)):
                    atoms.append(((i+fx)*a + (j+fy)*xy, (j+fy)*ly, (k+fz)*c))
    L = ["Ti","",f"{len(atoms)} atoms","1 atom types","",
         f"0 {lx:.6f} xlo xhi", f"0 {ly:.6f} ylo yhi", f"0 {lz:.6f} zlo zhi",
         f"{xy:.6f} 0.0 0.0 xy xz yz", "", "Masses", "", "1 47.867", "",
         "Atoms # atomic", ""]
    for i, p in enumerate(atoms, 1):
        L.append(f"{i} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    (W/"c.data").write_text("\n".join(L)+"\n", encoding="utf-8")
    inp = f"""units metal
atom_style atomic
read_data c.data
pair_style nep
pair_coeff * * {MODEL} Ti
neighbor 2.0 bin
neigh_modify one 4000 delay 0 every 1 check yes
thermo 100
min_style cg
min_modify dmax 0.02 line quadratic
minimize 1e-12 1e-12 4000 40000
print "EPER $(pe)"
"""
    (W/"in.eos").write_text(inp, encoding="utf-8")
    t0 = time.time()
    subprocess.run([EXE, "-in", "in.eos", "-log", "log.eos", "-screen", "none"],
                   cwd=W, check=True, stdout=subprocess.DEVNULL)
    dt = time.time() - t0
    txt = (W/"log.eos").read_text(errors="replace")
    m = re.findall(r"EPER\s+(-?\d+\.\d+)", txt)
    return (float(m[-1])/len(atoms) if m else None), dt

grid = [(round(a,3), round(cr,3))
        for a in np.arange(2.90, 3.03, 0.02)
        for cr in np.arange(1.55, 1.66, 0.025)]
out = {}
print(f"[Ti] {len(grid)} points, exact hexagonal cell")
best = None
for a, cr in grid:
    e, dt = energy(a, cr)
    if e is None:
        print(f"  a={a:.3f} c/a={cr:.3f}  FAILED"); continue
    out[f"{a}_{cr}"] = e
    if best is None or e < best[0]:
        best = (e, a, cr)
    print(f"  a={a:.3f} c/a={cr:.3f}  E={e:.4f}  ({dt:.0f}s)")
print(f"[MIN] a={best[1]:.3f} c/a={best[2]:.3f} E={best[0]:.4f}")
(W/"ti_grid.json").write_text(json.dumps({"grid": out,
    "min": {"a": best[1], "covera": best[2], "E": best[0]}}, indent=2))
