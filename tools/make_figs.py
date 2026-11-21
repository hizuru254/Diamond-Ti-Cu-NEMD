#!/usr/bin/env python3
"""Figures for the deliverable: temperature profiles and TBR ranking."""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/cc/si_ti_cu_nemd")
OUT = ROOT/"docs/figures"; OUT.mkdir(parents=True, exist_ok=True)

def profile(run, seed="seed_1001", nframes=100):
    rows, ts = [], None
    f = ROOT/"workspace"/run/seed/"profile.txt"
    for line in f.read_text().splitlines():
        s = line.split()
        if len(s) == 3 and all(t.isdigit() for t in s):
            ts = int(s[0]); continue
        if len(s) == 4 and s[0].isdigit() and ts is not None:
            try: rows.append((ts, float(s[1]), int(s[2]), float(s[3])))
            except ValueError: pass
    a = np.array(rows)
    late = np.unique(a[:,0])[-nframes:]
    a = a[np.isin(a[:,0], late) & (a[:,2] > 0)]
    z = np.unique(a[:,1])
    T = np.array([a[a[:,1]==zz, 3].mean() for zz in z])
    return z, T

fig, ax = plt.subplots(figsize=(6.4, 4.2))
for run, label in (("run_v005", "flat"), ("run_v003", "alloy 6 A"),
                   ("run_v002", "interlock H=6")):
    try:
        z, T = profile(run)
    except FileNotFoundError:
        continue
    m = (z > -55) & (z < 55)
    ax.plot(z[m], T[m], ".", ms=2, label=f"{label} ({run})")
ax.axvline(0, color="k", lw=0.6, ls="--")
ax.set_xlabel("z (A)"); ax.set_ylabel("T (K)")
ax.set_title("Temperature profiles (last 100 ps)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(OUT/"fig_profiles.png", dpi=160)
print("saved", OUT/"fig_profiles.png")

# ranking bar chart
data = json.loads((ROOT/"reports/pattern_ranking.json").read_text())
data = sorted([d for d in data if d.get("tbr")], key=lambda d: d["tbr"])
names = {"run_v005": "flat (control)", "run_v001": "flat (3 seeds)",
         "run_v003": "alloy 6 A", "run_v006": "interlock H=3",
         "run_v002": "interlock H=6", "run_v004": "alloy 10 A",
         "run_v007": "flat (repeat)"}
fig, ax = plt.subplots(figsize=(6.4, 3.8))
labels = [names.get(d["run"], d["run"]) for d in data]
vals = [d["tbr"]*1e10 for d in data]
cols = ["#2a9d8f" if v < 8 else "#e76f51" for v in vals]
ax.barh(range(len(vals)), vals, color=cols)
ax.set_yticks(range(len(vals))); ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("interface TBR (1e-10 m2 K / W)")
ax.set_title("Interface resistance: flat vs patterned")
for i, v in enumerate(vals):
    ax.text(v+0.2, i, f"{v:.2f}", va="center", fontsize=8)
ax.grid(axis="x", alpha=0.3)
fig.tight_layout(); fig.savefig(OUT/"fig_ranking.png", dpi=160)
print("saved", OUT/"fig_ranking.png")
