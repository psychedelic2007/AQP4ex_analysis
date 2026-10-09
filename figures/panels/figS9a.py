"""Fig. S9a - alpha and beta at tested extension codons in each robustness run."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

RUNS = [("MUSCLE", "aqp4_joint_muscle_extension"), ("PRANK", "aqp4_joint_prank_extension"),
        ("MAFFT dedup", "aqp4_joint_dedup_extension")] + \
       [(lab, f"aqp4_depth_{k}_joint_extension") for k, lab in C.DEPTH_LABEL.items()]
pos = list(range(325, 353))
fig, axs = plt.subplots(len(RUNS), 1, figsize=(7.1, 1.05 * len(RUNS)), sharex=True)
for ax, (lab, tab) in zip(axs, RUNS):
    t = C.site_table(tab)
    t = t[t.tested & t.human_residue.astype(str).str.isdigit()]
    xp = t.human_residue.astype(int)
    for xi, (_, r) in zip(xp, t.iterrows()):
        ax.plot([xi, xi], [max(r.beta, 1e-3), r.alpha], color=C.GRID, lw=1.2, zorder=1)
    ax.scatter(xp, t.alpha, s=20, facecolor="white", edgecolor=C.INK2, lw=0.6, zorder=2)
    ax.scatter(xp, t.beta.clip(lower=1e-3), s=20, color=C.ORANGE, lw=0, zorder=3)
    for xi, (_, r) in zip(xp, t.iterrows()):
        if r.purifying or r.diversifying:
            col = C.RED if r.diversifying else C.BLUE
            ax.scatter([xi], [6], marker="s", s=20, color=col if (r.bonf or r.diversifying) else "white",
                       edgecolor=col, lw=0.6, clip_on=False)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 8)
    ax.axhline(1, color=C.INK2, lw=0.4, ls="--")
    ax.text(1.005, 0.5, f"{lab}\n{int(t.n_taxa.iloc[0])} sp.", transform=ax.transAxes, fontsize=5.8,
            va="center")
    ax.tick_params(labelsize=10)
axs[-1].set_xticks(pos)
axs[-1].set_xticklabels(pos, rotation=90, fontsize=10)
axs[-1].set_xlabel("Human residue")
fig.supylabel("Rate (open α, filled β)", fontsize=6.5)
fig.tight_layout(h_pad=0.3)
C.save(fig, "FigS9", "a")
