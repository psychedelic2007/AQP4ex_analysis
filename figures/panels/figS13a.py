"""Fig. S13a - AQP5 restricted to the 124 selection-set species carrying UGA-C:
alpha and beta at tested downstream-frame columns (joint fit joint_ugac)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f = C.fel_sites("aqp5", "joint_ugac")
e = f[(f.module == "extension") & f.tested].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(3.2, 2.0))
x = np.arange(len(e))
for xi, (_, r) in zip(x, e.iterrows()):
    ax.plot([xi, xi], [r.beta, r.alpha], color=C.GRID, lw=1.4, zorder=1)
ax.scatter(x, e.alpha, s=14, facecolor="white", edgecolor=C.INK2, lw=0.7, zorder=2, label="α")
ax.scatter(x, e.beta, s=14, color=C.ORANGE, lw=0, zorder=3, label="β")
ax.set_yscale("log")
ax.axhline(1, color=C.INK2, lw=0.4, ls="--")
ax.set_xticks(x)
ax.set_xticklabels([f"col {c - 265}" for c in e.column], rotation=90, fontsize=5.3)
ax.set_ylabel("Rate")
ax.set_title(f"AQP5, UGA-C species (n = {int(f.n_taxa.iloc[0])}): {int(e.purifying.sum())}/{len(e)} "
             f"purifying\njoint global dN/dS {C.global_omega('aqp5', 'joint_ugac'):.3f}",
             fontweight="normal", fontsize=6)
ax.legend(fontsize=5.5, loc="lower right")
C.save(fig, "FigS13", "a")
