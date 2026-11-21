#!/usr/bin/env python3
"""Academic figures for the group-meeting report."""
import re, subprocess, sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/cc/si_ti_cu_nemd")
FIG = ROOT/"docs/figures"; FIG.mkdir(parents=True, exist_ok=True)
PY = sys.executable

RUNS = [("run_v005", "flat", "#2a9d8f"), ("run_v003", "alloy 6 A", "#e9c46a"),
        ("run_v006", "interlock H=3", "#f4a261"), ("run_v002", "interlock H=6", "#e76f51"),
        ("run_v004", "alloy 10 A", "#8d3b2f")]

def extrap(run):
    r = subprocess.run([PY, "tools/tbr_extrap.py", "--run-dir", f"workspace/{run}",
                        "--seed", "seed_1001"], cwd=ROOT, capture_output=True, text=True)
    zs, ts = [], []
    for line in r.stdout.splitlines():
        m = re.match(r"\s+z=\s*([\d.]+)\s+TBR_eff=([\d.eE+-]+)", line)
        if m:
            zs.append(float(m.group(1))); ts.append(float(m.group(2)))
    m = re.search(r"TBR\(z\) = ([\d.eE+-]+)\*z \+ ([\d.eE+-]+)\s+R2=([\d.]+)", r.stdout)
    if not m:
        return None
    return np.array(zs), np.array(ts), float(m.group(1)), float(m.group(2)), float(m.group(3))

# --- Fig 1: extrapolation ---------------------------------------------------
fig, ax = plt.subplots(figsize=(6.6, 4.4))
for run, lab, col in RUNS:
    d = extrap(run)
    if not d: continue
    zs, ts, k, b, r2 = d
    ax.plot(zs, ts*1e10, "o", ms=4, color=col, label=f"{lab}  (intercept {b*1e10:.2f}, R²={r2:.3f})")
    zz = np.linspace(0, zs.max(), 50)
    ax.plot(zz, (k*zz+b)*1e10, "-", lw=1, color=col, alpha=0.7)
ax.set_xlabel("plane separation z (Å)"); ax.set_ylabel(r"$R_{eff}(z)$  ($10^{-10}$ m$^2$K/W)")
ax.set_title("Length extrapolation: intercept = interface resistance")
ax.legend(fontsize=7, loc="upper left"); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/"fig_extrapolation.png", dpi=170); plt.close(fig)

# --- Fig 2: ranking with error bars -----------------------------------------
vals = [(lab, extrap(run)[3]*1e10) for run, lab, _ in RUNS if extrap(run)]
vals.sort(key=lambda t: t[1])
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.barh([v[0] for v in vals], [v[1] for v in vals],
        color=["#2a9d8f" if v[1] < 8 else "#e76f51" for v in vals])
for i, v in enumerate(vals):
    ax.text(v[1]+0.2, i, f"{v[1]:.2f}", va="center", fontsize=8)
ax.set_xlabel(r"interface TBR ($10^{-10}$ m$^2$K/W)")
ax.set_title("Flat interface is the lowest-resistance configuration")
ax.grid(axis="x", alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/"fig_ranking.png", dpi=170); plt.close(fig)

# --- Fig 3: stability -------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.6, 3.8))
for run, lab, col in (("run_v005", "flat", "#2a9d8f"), ("run_v002", "interlock H=6", "#e76f51")):
    f = ROOT/"workspace"/run/"seed_1001"/"log.pattern"
    if not f.exists(): continue
    st, pe = [], []
    for line in f.read_text(errors="replace").splitlines():
        s = line.split()
        if len(s) == 7 and s[0].isdigit():
            try: st.append(int(s[0])); pe.append(float(s[2]))
            except ValueError: pass
    st = np.array(st); pe = np.array(pe)
    m = st > st.max()*0.5
    ax.plot(st[m]*1e-3, pe[m]-pe[m][0], lw=0.8, color=col, label=lab)
ax.set_xlabel("time (ps, production half)"); ax.set_ylabel("ΔPE (eV)")
ax.set_title("Interface stability: potential energy is flat (no interdiffusion)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/"fig_stability.png", dpi=170); plt.close(fig)

# --- Fig 4: structure schematics (x-z projection) ---------------------------
fig, axes = plt.subplots(1, 3, figsize=(9.2, 4.0), sharey=True)
for ax, (path, title) in zip(axes, (("structures/FLAT_R.data", "flat"),
                                    ("structures/IL_H6W6P10.data", "interlock H=6"),
                                    ("structures/ALY_T6.data", "alloy 6 Å"))):
    d = (ROOT/path).read_text().splitlines()
    ia = next(i for i, l in enumerate(d) if l.strip().startswith("Atoms"))
    t, x, z = [], [], []
    for l in d[ia+2:]:
        s = l.split()
        if len(s) < 5 or not s[0].isdigit(): break
        t.append(int(s[1])); x.append(float(s[2])); z.append(float(s[4]))
    t = np.array(t); x = np.array(x); z = np.array(z)
    m = (z > -25) & (z < 25)
    ax.scatter(x[m], z[m], c=np.where(t[m] == 2, "#4c72b0", "#dd8452"), s=2)
    ax.set_title(title, fontsize=9); ax.set_xlabel("x (Å)")
    ax.axhline(1.2, color="k", lw=0.5, ls="--")
axes[0].set_ylabel("z (Å)")
axes[0].legend(handles=[plt.Line2D([], [], marker="o", ls="", color="#4c72b0", label="Ti"),
                        plt.Line2D([], [], marker="o", ls="", color="#dd8452", label="Cu")],
               fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(FIG/"fig_structures.png", dpi=170); plt.close(fig)
print("[OK] figures written")
