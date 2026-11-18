#!/usr/bin/env python3
"""Parallel hcp Ti grid (fixed cell geometry): per-cell y extent, correct basis."""
import json, re, subprocess, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

BASE = Path("/tmp/ti_grid3"); BASE.mkdir(parents=True, exist_ok=True)
EXE = "/home/cc/lammps/src/lmp_mpi"
MODEL = "/home/cc/si_ti_cu_nemd/potentials/CuTiC_NEP.txt"
NX, NY, NZ = 3, 3, 2

def build(a, cr, d):
    c = a * cr
    cy = a * np.sqrt(3) / 2.0          # y extent of ONE hexagonal cell
    lx, ly, lz = NX*a, NY*cy, NZ*c
    xy = -a/2
    atoms = []
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                for fx, fy, fz in ((0,0,0), (1/3,2/3,.5)):
                    atoms.append(((i+fx)*a + (j+fy)*xy, (j+fy)*cy, (k+fz)*c))
    L = ["Ti","",f"{len(atoms)} atoms","1 atom types","",
         f"0 {lx:.6f} xlo xhi", f"0 {ly:.6f} ylo yhi", f"0 {lz:.6f} zlo zhi",
         f"{xy:.6f} 0.0 0.0 xy xz yz", "", "Masses", "", "1 47.867", "",
         "Atoms # atomic", ""]
    for i, p in enumerate(atoms, 1):
        L.append(f"{i} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    (d/"c.data").write_text("\n".join(L)+"\n", encoding="utf-8")
    return len(atoms)

def run_point(args):
    a, cr = args
    d = BASE / f"a{a:.3f}_c{cr:.3f}"; d.mkdir(exist_ok=True)
    n = build(a, cr, d)
    inp = f"""units metal
atom_style atomic
read_data c.data
pair_style nep
pair_coeff * * {MODEL} Ti
neighbor 2.0 bin
neigh_modify one 4000 delay 0 every 1 check yes
min_style cg
min_modify dmax 0.01 line quadratic
minimize 1e-12 1e-12 2000 20000
print "EPER $(pe)"
"""
    (d/"in.eos").write_text(inp, encoding="utf-8")
    t0 = time.time()
    try:
        subprocess.run([EXE,"-in","in.eos","-log","log.eos","-screen","none"],
                       cwd=d, timeout=120, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.TimeoutExpired:
        return (a, cr, None, time.time()-t0, "timeout")
    txt = (d/"log.eos").read_text(errors="replace")
    m = re.findall(r"EPER\s+(-?\d+\.\d+)", txt)
    return (a, cr, (float(m[-1])/n if m else None), time.time()-t0, "ok" if m else "noe")

grid = [(round(a,3), round(cr,3))
        for a in np.arange(2.90, 3.03, 0.02)
        for cr in np.arange(1.55, 1.66, 0.025)]
print(f"[Ti] {len(grid)} points, 24 workers", flush=True)
t0 = time.time(); rows = []
with ThreadPoolExecutor(max_workers=24) as ex:
    for r in ex.map(run_point, grid):
        rows.append(r); a, cr, e, dt, st = r
        print(f"  a={a:.3f} c/a={cr:.3f}  E={('%.4f'%e) if e is not None else '--':>8s} ({dt:.0f}s {st})", flush=True)
ok = [r for r in rows if r[2] is not None]; ok.sort(key=lambda r: r[2])
print(f"[DONE] {len(ok)}/{len(rows)} in {time.time()-t0:.0f}s")
for r in ok[:8]:
    print(f"  low: a={r[0]:.3f} c/a={r[1]:.3f} E={r[2]:.4f}")
json.dump({"rows": [[r[0],r[1],r[2],r[4]] for r in rows]}, open(BASE/"grid.json","w"), indent=2)
