import subprocess, re, json
from pathlib import Path
import numpy as np
W = Path("/tmp/tichk"); W.mkdir(parents=True, exist_ok=True)
EXE="/home/cc/lammps/src/lmp_mpi"; MODEL="/home/cc/si_ti_cu_nemd/potentials/nep89.txt"

def build(kind, a, covera=None, nx=3, ny=3, nz=2):
    if kind == "hcp":
        c = a*covera; lx, ly, lz = a*nx, a*np.sqrt(3)*ny, c*nz
        basis = [(0,0,0), (2/3, 1/3, 0.5)]
        pts = [((i+bx)*a, (j+by)*a*np.sqrt(3), (k+bz)*c)
               for i in range(nx) for j in range(ny) for k in range(nz) for bx,by,bz in basis]
    elif kind == "fcc":
        lx = ly = lz = a*nx
        basis = [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)]
        pts = [((i+bx)*a,(j+by)*a,(k+bz)*a) for i in range(nx) for j in range(ny) for k in range(nz) for bx,by,bz in basis]
    else:  # bcc
        lx = ly = lz = a*nx
        basis=[(0,0,0),(.5,.5,.5)]
        pts = [((i+bx)*a,(j+by)*a,(k+bz)*a) for i in range(nx) for j in range(ny) for k in range(nz) for bx,by,bz in basis]
    lines=[kind,"",f"{len(pts)} atoms","1 atom types","",f"0 {lx:.6f} xlo xhi",f"0 {ly:.6f} ylo yhi",f"0 {lz:.6f} zlo zhi","","Masses","","1 47.867","","Atoms # atomic",""]
    for i,p in enumerate(pts,1): lines.append(f"{i} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    (W/"c.data").write_text("\n".join(lines)+"\n", encoding="utf-8")
    return len(pts), (lx,ly,lz)

def relax(kind, a, covera=None):
    n, box = build(kind, a, covera)
    inp=f"""units metal
atom_style atomic
boundary p p p
read_data c.data
pair_style nep
pair_coeff * * {MODEL} Ti
neighbor 2.0 bin
thermo_style custom step pe pxx pyy pzz
fix br all box/relax aniso 0.0 vmax 0.0005
neigh_modify one 4000 delay 0 every 1 check yes
minimize 1e-12 1e-12 20000 200000
variable pe equal pe
print "RESULT $(pe) $(lx) $(ly) $(lz)"
"""
    (W/"in.r").write_text(inp, encoding="utf-8")
    r=subprocess.run([EXE,"-in","in.r","-log","log.r","-screen","none"],cwd=W,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    txt=(W/"log.r").read_text(errors="replace")
    m=re.findall(r"RESULT\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)",txt)
    if not m: return None, txt[-300:]
    pe,lx,ly,lz=[float(x) for x in m[-1]]
    return (pe/n, n, (lx,ly,lz)), ""

print("=== relax from experimental-ish cells (box/relax tri) ===")
for kind,a,cr in (("hcp",2.95,1.587),("fcc",4.10,None),("bcc",3.30,None)):
    r,err = relax(kind,a,cr)
    if r is None:
        print(f"{kind}: FAILED {err}"); continue
    e,n,box = r
    if kind=="hcp":
        lx,ly,lz=box; nx=3; ny=3; nz=2
        aa=lx/nx; cc=lz/nz
        print(f"hcp: E={e:.4f} eV/atom  a={aa:.4f} c={cc:.4f} c/a={cc/aa:.3f}")
    else:
        lx,ly,lz=box; aa=lx/3
        print(f"{kind}: E={e:.4f} eV/atom  a={aa:.4f}")
