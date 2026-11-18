#!/usr/bin/env python3
"""Relax hcp Ti with the dedicated NEP from the experimental cell, larger supercell."""
import re, subprocess
from pathlib import Path
import numpy as np
W = Path("/tmp/ti2"); W.mkdir(parents=True, exist_ok=True)
EXE = "/home/cc/lammps/src/lmp_mpi"; MODEL = "/home/cc/si_ti_cu_nemd/potentials/CuTiC_NEP.txt"

def build(a, covera, n=4, nz=3):
    c = a*covera
    lx, ly, lz = a*n, a*np.sqrt(3)*n, c*nz
    basis = [(0,0,0), (2/3,1/3,.5)]
    pts = [((i+bx)*a, (j+by)*a*np.sqrt(3), (k+bz)*c)
           for i in range(n) for j in range(n) for k in range(nz) for bx,by,bz in basis]
    L = ["Ti","",f"{len(pts)} atoms","1 atom types","",f"0 {lx:.6f} xlo xhi",
         f"0 {ly:.6f} ylo yhi", f"0 {lz:.6f} zlo zhi","","Masses","","1 47.867","",
         "Atoms # atomic",""]
    for i,p in enumerate(pts,1): L.append(f"{i} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    (W/"c.data").write_text("\n".join(L)+"\n", encoding="utf-8")
    return len(pts), (lx,ly,lz)

for (a0, cr0) in ((2.951,1.587), (2.951,1.500), (2.90,1.63)):
    n, box = build(a0, cr0)
    inp = f"""units metal
atom_style atomic
boundary p p p
read_data c.data
pair_style nep
pair_coeff * * {MODEL} Ti
neighbor 2.0 bin
neigh_modify one 4000
thermo_style custom step pe pxx pyy pzz
fix br all box/relax aniso 0.0 vmax 0.0002
minimize 1e-12 1e-12 50000 500000
print "BOX $(lx) $(ly) $(lz) $(pe)"
"""
    (W/"in.r").write_text(inp, encoding="utf-8")
    subprocess.run([EXE,"-in","in.r","-log","log.r","-screen","none"], cwd=W,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t = (W/"log.r").read_text(errors="replace")
    m = re.findall(r"BOX\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", t)
    if not m:
        print(f"start a={a0} c/a={cr0}: FAILED"); continue
    lx,ly,lz,pe = [float(x) for x in m[-1]]
    aa = lx/4; cc = lz/3
    print(f"start a={a0:.3f} c/a={cr0:.3f}  ->  a={aa:.4f} c={cc:.4f} c/a={cc/aa:.3f} "
          f"E={pe/n:.4f} eV/atom")
