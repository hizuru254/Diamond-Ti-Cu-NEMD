#!/usr/bin/env python3
"""Credibility figure: running mean convergence + spread of the seven measurements."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/cc/si_ti_cu_nemd"); FIG = ROOT/"docs/figures"
sys.path.insert(0, str(ROOT/"tools"))
import geomlib
from tbr_plane import parse_blocks, profile, dT_planes

def series(run, nblocks=15, z=30.0):
    r = ROOT/"workspace"/run
    geom = geomlib.load_json(r/"geom.json"); params = geomlib.load_json(r/"params.snapshot.json")
    blocks = parse_blocks(r/"seed_1001"/"profile.txt")
    zc = float(geom["interfaces"][0]["z"])
    area = float(geom["box"]["lx"]) * float(geom["box"]["ly"])
    q = (float(params["flux_ev_per_ps"])*float(params["units"]["ev_per_ps_to_W"])
         / (area*float(params["units"]["ang2_to_m2"])))
    idx = np.linspace(0, len(blocks), nblocks+1).astype(int)
    v = []
    for a, b in zip(idx[:-1], idx[1:]):
        if b <= a: continue
        zb, Tb = profile(blocks[a:b]); v.append(dT_planes(zb, Tb, z, zc)/q)
    return np.array(v)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 3.9), gridspec_kw={"width_ratios": [1.35, 1]})

v = series("run_v005")
cum = np.cumsum(v)/np.arange(1, len(v)+1)
ax1.plot(np.arange(1, len(v)+1), v*1e10, "o", ms=3.5, color="#9ec9b6", label="block value")
ax1.plot(np.arange(1, len(v)+1), cum*1e10, "-", lw=2, color="#1F5C43", label="running mean")
ax1.axhline(cum[-1]*1e10, ls="--", lw=0.8, color="#1F5C43", alpha=0.6)
ax1.set_xlabel("block (of the 150 ps production run)")
ax1.set_ylabel(r"$R_{eff}(\pm30\AA)$  ($10^{-10}$ m$^2$K/W)")
ax1.set_title("Running mean converges: %.1f $\\pm$ %.1f" % (cum[-1]*1e10, v.std(ddof=1)/np.sqrt(len(v))*1e10),
               fontsize=11)
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

vals = np.array([7.292, 7.464, 7.229, 6.931, 6.931, 7.642, 6.290])
lab = ["s1", "s2", "s3", "rebuild", "repeat", "2 eV/ps", "3 eV/ps"]
ax2.axhspan(vals.mean()-vals.std(ddof=1), vals.mean()+vals.std(ddof=1),
            color="#E8F4EC", zorder=0, label="±1σ")
ax2.plot(vals, "o", ms=7, color="#1F5C43")
for i, (l, x) in enumerate(zip(lab, vals)):
    ax2.annotate(l, (i, x), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=7.5)
ax2.axhline(vals.mean(), ls="--", lw=1, color="#1F5C43")
ax2.set_xticks([]); ax2.set_ylabel(r"$R_{int}$  ($10^{-10}$ m$^2$K/W)")
ax2.set_title("Seven independent measurements: %.1f $\\pm$ %.1f" %
              (vals.mean(), vals.std(ddof=1)), fontsize=11)
ax2.set_xlim(-1, len(vals)); ax2.legend(fontsize=8, loc="lower left"); ax2.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(FIG/"fig_credibility.png", dpi=170)
print("[OK] fig_credibility.png")
