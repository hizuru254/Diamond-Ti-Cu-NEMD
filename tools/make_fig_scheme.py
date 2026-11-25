#!/usr/bin/env python3
"""Pattern schematic: flat vs interlocked, with H/W/P labelled."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

FIG = Path("/home/cc/si_ti_cu_nemd/docs/figures")
TI = "#4C72B0"; CU = "#DD8452"

fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4))

# left: flat
ax = axes[0]
ax.add_patch(Rectangle((0, 0), 10, 1.6, color=TI))
ax.add_patch(Rectangle((0, 1.6), 10, 2.2, color=CU))
ax.axhline(1.6, color="k", lw=0.8, ls="--")
ax.text(5, 2.7, "Cu", ha="center", va="center", color="w", fontsize=13, weight="bold")
ax.text(5, 0.8, "Ti", ha="center", va="center", color="w", fontsize=13, weight="bold")
ax.set_title("flat reference", fontsize=12)
ax.set_xlim(0, 10); ax.set_ylim(0, 3.8); ax.axis("off")

# right: interlocked
ax = axes[1]
H, W, P = 0.55, 1.2, 2.5
ax.add_patch(Rectangle((0, 0), 10, 1.6, color=TI))
xs = 0.0
while xs < 10:
    ax.add_patch(Rectangle((xs, 1.6), W, H, color=TI))      # ridge
    xs += P
ax.add_patch(Rectangle((0, 1.6 + H), 10, 2.2 - H, color=CU))
# fill trenches with Cu
xs = W
while xs < 10:
    ax.add_patch(Rectangle((xs, 1.6), min(P - W, 10 - xs), H, color=CU))
    xs += P
# labels
ax.annotate("", xy=(0.0, 1.6), xytext=(0.0, 1.6 + H),
            arrowprops=dict(arrowstyle="<->", lw=1))
ax.text(-0.35, 1.6 + H/2, "H", ha="right", va="center", fontsize=11, color="#1F5C43", weight="bold")
ax.annotate("", xy=(0.0, 2.45), xytext=(W, 2.45),
            arrowprops=dict(arrowstyle="<->", lw=1))
ax.text(W/2, 2.55, "W", ha="center", va="bottom", fontsize=11, color="#1F5C43", weight="bold")
ax.annotate("", xy=(0.0, 3.15), xytext=(P, 3.15),
            arrowprops=dict(arrowstyle="<->", lw=1, color="#666"))
ax.text(P/2, 3.22, "P", ha="center", va="bottom", fontsize=11, color="#666", weight="bold")
ax.set_title("interlocked: Ti ridges, Cu fills the trenches", fontsize=11)
ax.set_xlim(-0.8, 10); ax.set_ylim(0, 3.8); ax.axis("off")

fig.tight_layout(); fig.savefig(FIG/"fig_pattern_schematic.png", dpi=180)
print("[OK] fig_pattern_schematic.png")
