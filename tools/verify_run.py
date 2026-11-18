#!/usr/bin/env python3
"""Independent verification: profile, PE drift, interdiffusion (clean parser)."""
import sys
from pathlib import Path
import numpy as np

run = Path(sys.argv[1])
sd = sys.argv[2] if len(sys.argv) > 2 else "seed_1001"
seed = run / sd
print(f"[verify] {run.name}/{sd}")
rows = []
ts = None
for line in (seed/"profile.txt").read_text().splitlines():
    s = line.split()
    if len(s) == 3 and all(t.isdigit() for t in s):
        ts = int(s[0]); continue
    if len(s) == 4 and s[0].isdigit() and ts is not None:
        try:
            rows.append((ts, int(s[0]), float(s[1]), int(s[2]), float(s[3])))
        except ValueError:
            continue
a = np.array(rows, dtype=float)
allts = np.unique(a[:,0])
late = allts[-min(100, len(allts)):]
sel = a[np.isin(a[:,0], late)]
sel = sel[sel[:,3] > 0]
print(f"[profile] {len(allts)} frames; averaging last {len(late)}; {len(sel)} filled bins")
z, T = sel[:,2], sel[:,4]
for lo in range(-60, 61, 10):
    m = (z >= lo) & (z < lo+10)
    if m.sum():
        print(f"  z {lo:+4d}..{lo+10:+4d}  T={T[m].mean():7.2f} K  (n={m.sum()})")

log = (seed/"log.pattern").read_text().splitlines()
pe = []
for line in log:
    s = line.split()
    if len(s) == 7 and s[0].isdigit():
        try: pe.append((int(s[0]), float(s[1]), float(s[2])))
        except ValueError: pass
if pe:
    pe = np.array(pe)
    half = pe[len(pe)//2:]
    slope = np.polyfit(half[:,0], half[:,2], 1)[0] * 1000
    print(f"[energy] {len(pe)} thermo rows; production-half dPE/dt = {slope:+.4f} eV/ps")

d = (seed/"final.data").read_text().splitlines()
ia = next(i for i,l in enumerate(d) if l.strip().startswith("Atoms"))
ty, zz = [], []
for l in d[ia+2:]:
    s = l.split()
    if len(s) < 5 or not s[0].isdigit(): break
    ty.append(int(s[1])); zz.append(float(s[4]))
ty, zz = np.array(ty), np.array(zz)
ti, cu = zz[ty==2], zz[ty==3]
print(f"[structure] Ti z {ti.min():.1f}..{ti.max():.1f} | Cu z {cu.min():.1f}..{cu.max():.1f}")
print(f"[interdiffusion] Ti atoms with z > interface+5: {(ti > 6.2).sum()} "
      f"| Cu atoms with z < interface-5: {(cu < -3.8).sum()}")
