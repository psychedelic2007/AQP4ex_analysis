"""Fig. 4d - alpha vs beta at every tested extension/downstream-frame codon, per gene.
Filled = purifying (q <= 0.20). AQP4 purifying sites are repeated in grey in each paralogue
facet for reference; inset text gives median alpha in the downstream frame vs channel."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

fig, axs = plt.subplots(2, 3, figsize=(5.4, 3.7), sharex=True, sharey=True)
ref = None
for ax, g in zip(axs.flat, C.GENES):
    f = C.fel_sites(g, "joint_mafft")
    e = f[(f.module == "extension") & f.tested]
    c = f[(f.module == "channel") & f.tested]
    if g == "aqp4":
        ref = e[e.purifying]
    else:
        ax.scatter(ref.alpha, ref.beta, s=8, color=C.NEUTRAL, lw=0, alpha=0.7,
                   label="AQP4 purifying sites")
    col = C.BLUE if g == "aqp4" else C.ORANGE
    ax.scatter(e[~e.purifying].alpha, e[~e.purifying].beta, s=12, facecolor="white",
               edgecolor=col, lw=0.7, label="not called")
    ax.scatter(e[e.purifying].alpha, e[e.purifying].beta, s=14, color=col,
               edgecolor="white", lw=0.3, label="purifying, q ≤ 0.20")
    ax.plot([0.005, 20], [0.005, 20], color=C.INK2, lw=0.5, ls="--")
    ax.axvline(1, color=C.GRID, lw=0.6, zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 12)
    ax.set_ylim(0.008, 12)
    ax.set_title(f"{C.GENE_LABEL[g]}  ({int(e.purifying.sum())}/{len(e)} purifying)",
                 fontweight="normal", fontsize=6.5)
    ax.text(0.97, 0.04, f"median α\next {e.alpha.median():.2f}\nchannel {c.alpha.median():.2f}",
            transform=ax.transAxes, ha="right", fontsize=5.3, color=C.INK2)
for ax in axs[1]:
    ax.set_xlabel("α (synonymous)")
for ax in axs[:, 0]:
    ax.set_ylabel("β (non-synonymous)")
axs[0, 1].legend(loc="upper left", fontsize=5, handletextpad=0.1, borderaxespad=0.2)
fig.tight_layout(h_pad=0.8, w_pad=0.5)
C.save(fig, "Fig4", "d")
