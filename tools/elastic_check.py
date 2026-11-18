import subprocess, re, json
from pathlib import Path
import numpy as np
W=Path("/tmp/elast"); W.mkdir(parents=True, exist_ok=True)
EXE="/home/cc/lammps/src/lmp_mpi"; MODEL="/home/cc/si_ti_cu_nemd/potentials/nep89.txt"
MASS={"Si":28.0855,"Cu":63.546,"Ti":47.867}
BAS={"Si":[(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0),(.25,.25,.25),(.25,.75,.75),(.75,.25,.75),(.75,.75,.25)],
     "Cu":[(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)]}

def cell(el, a, c=None, n=2):
    if el=="Ti":
        c=c or a*1.414
        basis=[(0,0,0),(2/3,1/3,.5)]
        lx,ly,lz=a*n, a*np.sqrt(3)*n, c*n
        pts=[((i+bx)*a,(j+by)*a*np.sqrt(3),(k+bz)*c) for i in range(n) for j in range(n) for k in range(n) for bx,by,bz in basis]
    else:
        basis=BAS[el]; lx=ly=lz=a*n
        pts=[((i+bx)*a,(j+by)*a,(k+bz)*a) for i in range(n) for j in range(n) for k in range(n) for bx,by,bz in basis]
    L=[el,"",f"{len(pts)} atoms","1 atom types","",f"0 {lx:.6f} xlo xhi",f"0 {ly:.6f} ylo yhi",f"0 {lz:.6f} zlo zhi","","Masses","",f"1 {MASS[el]}","","Atoms # atomic",""]
    for i,p in enumerate(pts,1): L.append(f"{i} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    (W/"c.data").write_text("\n".join(L)+"\n",encoding="utf-8")
    return len(pts), lx*ly*lz

def E(el,a,c=None,minimize=True):
    n,V=cell(el,a,c)
    cmd="minimize 1e-12 1e-12 20000 200000" if minimize else "run 0"
    inp=f"""units metal
atom_style atomic
boundary p p p
read_data c.data
pair_style nep
pair_coeff * * {MODEL} {el}
neighbor 2.0 bin
neigh_modify one 4000
thermo_style custom step pe
{cmd}
print "E $(pe)"
"""
    (W/"in.e").write_text(inp,encoding="utf-8")
    subprocess.run([EXE,"-in","in.e","-log","log.e","-screen","none"],cwd=W,check=True,stdout=subprocess.DEVNULL)
    t=(W/"log.e").read_text(errors="replace")
    m=re.findall(r"E\s+(-?\d+\.\d+)",t)
    return float(m[-1]), V

def bulk(el,a0,c0=None):
    Vs=[];Es=[]
    for s in (0.98,0.99,1.0,1.01,1.02):
        if el=="Ti":
            a=a0*s**(1/3); c=c0*s**(1/3)
        else:
            a=a0*s**(1/3); c=None
        e,V=E(el,a,c); Vs.append(V); Es.append(e)
    Vs=np.array(Vs); Es=np.array(Es)
    co=np.polyfit(Vs,Es,2)
    B=2*co[0]*Vs.min()   # d2E/dV2 * V
    ev=1.602176634e-19; A=1e-30
    return 2*co[0]*np.mean(Vs)*ev/A/1e9, np.polyfit(Vs,Es,2)[1]

print("=== bulk modulus from E(V), 5 volumes +-2% ===")
data={}
for el,a0,c0 in (("Si",5.400,None),("Cu",3.560,None),("Ti",2.859,4.043)):
    B,_=bulk(el,a0,c0)
    data[el]={"B_GPa":round(B,1)}
    print(f"{el}: B = {B:6.1f} GPa")
json.dump(data,open("/home/cc/si_ti_cu_nemd/potentials/elastic.json","w"),indent=2)
