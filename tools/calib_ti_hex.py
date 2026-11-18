#!/usr/bin/env python3
"""Exact hcp Ti cell (triclinic hexagonal), grid + convexity check, both NEP models."""
import re, subprocess
from pathlib import Path
import numpy as np
W = Path("/tmp/ti5"); W.mkdir(parents=True, exist_ok=True)
EXE = "/home/cc/lammps/src/lmp_mpi"

def energy(model, el, a, cr, nx=3, ny=3, nz=2):
    c = a*cr
    # triclinic hexagonal: a1=(a,0,0) a2=(-a/2, sqrt3*a/2, 0) c=(0,0,c)
    lx = nx*a; ly = ny*a*np.sqrt(3)/2; lz = nz*c; xy = -a/2
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for (fx, fy, fz) in ((0,0,0), (1/3, 2/3, .5)):
                    atoms.append(((i+fx)*a + (j+fy)*xy, (j+fy)*ly, (k+fz)*c))
    L = [el, "", f"{len(atoms)} atoms", "1 atom types", "",
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
pair_coeff * * {model} {el}
neighbor 2.0 bin
neigh_modify one 4000
thermo_style custom step pe
minimize 1e-12 1e-12 5000 50000
print "EPER $(pe)"
"""
    (W/"in.e").write_text(inp, encoding="utf-8")
    subprocess.run([EXE,"-in","in.e","-log","log.e","-screen","none"], cwd=W,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t = (W/"log.e").read_text(errors="replace")
    m = re.findall(r"EPER\s+(-?\d+\.\d+)", t)
    return float(m[-1])/len(atoms) if m else None

for tag, model in (("CuTiC_NEP", "/home/cc/si_ti_cu_nemd/potentials/CuTiC_NEP.txt"),
                   ("NEP89",     "/home/cc/si_ti_cu_nemd/potentials/nep89.txt")):
    print(f"=== {tag}  (exp: a=2.951 c/a=1.587 E0=-4.85) ===")
    print("        " + "  ".join(f"c/a={cr}" for cr in (1.55,1.587,1.62)))
    best=None
    for a in (2.92, 2.95, 2.98, 3.01):
        row=[]
        for cr in (1.55, 1.587, 1.62):
            e = energy(model, "Ti", a, cr); row.append(e)
            if e is not None and (best is None or e < best[0]): best = (e, a, cr)
        print(f"  a={a:.2f} " + "  ".join(f"{('%.3f'%v) if v else 'FAIL':>8s}" for v in row))
    if best:
        print(f"  -> min a={best[1]:.2f} c/a={best[2]:.3f} E={best[0]:.3f} eV/atom")
