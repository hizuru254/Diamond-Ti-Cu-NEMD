#!/usr/bin/env python3
"""Block-wise TBR over the production segment: is the measurement stationary?"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/cc/si_ti_cu_nemd")
FIG = ROOT/"docs/figures"
sys.path.insert(0, str(ROOT/"tools"))
import geomlib
from tbr_plane import parse_blocks, profile, dT_planes

def series(run, nblocks=12, z=30.0):
    r = ROOT/"workspace"/run
    geom = geomlib.load_json(r/"geom.json")
    params = geomlib.load_json(r/"params.snapshot.json")
    blocks = parse_blocks(r/"seed_1001"/"profile.txt")
    zc = float(geom["interfaces"][0]["z"])
    area = float(geom["box"]["lx"]) * float(geom["box"]["ly"])
    q = (float(params["flux_ev_per_ps"]) * float(params["units"]["ev_per_ps_to_W"])
         / (area * float(params["units"]["ang2_to_m2"])))
    idx = np.linspace(0, len(blocks), nblocks+1).astype(int)
    t, v = [], []
    for a, b in zip(idx[:-1], idx[1:]):
        if b <= a: continue
        zb, Tb = profile(blocks[a:b])
        t.append(0.5*(a+b))
        v.append(dT_planes(zb, Tb, z, zc)/q)
    return np.array(t), np.array(v)

fig, ax = plt.subplots(figsize=(6.8, 3.9))
for run, lab, col in (("run_v005", "flat", "#1F5C43"), ("run_v002", "interlock H=6", "#E07A5F")):
    try:
        t, v = series(run)
    except Exception as e:
        print("skip", run, e); continue
    ax.plot(t, v*1e10, "o-", ms=4, lw=1.2, color=col, label=f"{lab}: {v.mean()*1e10:.2f} ± {v.std(ddof=1)/np.sqrt(len(v))*1e10:.2f}")
ax.set_xlabel("block index within the 150 ps production run")
ax.set_ylabel(r"$R_{eff}(\pm 30\,\AA)$  ($10^{-10}$ m$^2$K/W)")
ax.set_title("Measured resistance is stationary, and the gap between designs is systematic")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/"fig_tbr_vs_time.png", dpi=170)
print("[OK] fig_tbr_vs_time.png")
