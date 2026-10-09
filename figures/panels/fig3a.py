"""Fig. 3a - Site-wise synonymous (alpha) and non-synonymous (beta) rates at the 28
tested extension codons (primary joint fit: MAFFT, 524 species)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from matplotlib.lines import Line2D  # noqa: E402
t = C.site_table("aqp4_joint_mafft_extension")
t = t[t.tested].copy()
t["pos"] = t.human_residue.astype(int)
x = np.arange(len(t))
SEG = [(325, 334, C.AQUA), (335, 344, C.NEUTRAL), (345, 352, C.BLUE)]

fig, ax = plt.subplots(figsize=(7.1, 2.3))
for a, b, col in SEG:
    i0, i1 = t.pos.tolist().index(a), t.pos.tolist().index(b)
    ax.axvspan(i0 - 0.5, i1 + 0.5, color=col, alpha=0.08, lw=0)
ax.axhline(1, color=C.INK2, lw=0.6, ls="--")
ax.text(len(t) - 0.4, 1.08, "gene-wide mean α = 1", ha="right", fontsize=5.5, color=C.INK2)
for xi, (_, r) in zip(x, t.iterrows()):
    ax.plot([xi, xi], [r.beta, r.alpha], color=C.GRID, lw=1.6, zorder=1)
ax.scatter(x, t.alpha, s=30, facecolor="white", edgecolor=C.INK2, lw=0.9, zorder=3)
ax.scatter(x, t.beta, s=30, color=C.ORANGE, edgecolor="white", lw=0.4, zorder=4)
ax.set_yscale("log")
ax.set_ylim(0.008, 4)
ax.set_ylabel("Substitution rate")
ax.set_xticks(x)
ax.set_xticklabels([f"{a}{p}" for a, p in zip(t.human_aa, t.pos)], rotation=90, fontsize=10)
ax.tick_params(axis="both", labelsize=10)
ax.set_xlim(-0.6, len(t) - 0.4)
# call track
for xi, (_, r) in zip(x, t.iterrows()):
    if r.purifying:
        ax.scatter([xi], [2.9], marker="s", s=30, color=C.BLUE if r.bonf else "white",
                   edgecolor=C.BLUE, lw=0.8, zorder=5, clip_on=False)
ax.text(-0.9, 2.9, "call", ha="right", va="center", fontsize=5.8)
h = [Line2D([], [], marker="o", ls="", mfc="white", mec=C.INK2, label="α (synonymous)"),
     Line2D([], [], marker="o", ls="", mfc=C.ORANGE, mec="white", label="β (non-synonymous)"),
     Line2D([], [], marker="s", ls="", mfc="white", mec=C.BLUE, label="purifying, q ≤ 0.20"),
     Line2D([], [], marker="s", ls="", mfc=C.BLUE, mec=C.BLUE, label="also Bonferroni (P ≤ 0.05/28)")]
ax.legend(handles=h, ncol=5, loc="lower left", bbox_to_anchor=(0, 1.08), fontsize=5.8,
          handletextpad=0.2, columnspacing=1.0)
C.recessive_grid(ax)
C.save(fig, "Fig3", "a")
