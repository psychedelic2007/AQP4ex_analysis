"""Fig. 4c - Module-specific global dN/dS (separate channel-only and extension-only fits)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

d = pd.DataFrame([dict(gene=g, ch=C.global_omega(g, "channel"), ex=C.global_omega(g, "ext_mafft"))
                  for g in C.GENES])
fig, ax = plt.subplots(figsize=(2.6, 2.0))
x = np.arange(len(d))
for i, r in d.iterrows():
    ax.plot([i, i], [r.ch, r.ex], color=C.GRID, lw=2, zorder=1)
ax.scatter(x, d.ch, s=22, color=C.INK2, label="channel", zorder=3, edgecolor="white", lw=0.4)
ax.scatter(x, d.ex, s=22, color=[C.BLUE if g == "aqp4" else C.ORANGE for g in d.gene],
           label="downstream frame", zorder=3, edgecolor="white", lw=0.4)
for i, r in d.iterrows():
    ax.text(i + 0.14, r.ex, f"{r.ex:.2f}", va="center", fontsize=5.5)
    ax.text(i + 0.14, r.ch, f"{r.ch:.2f}", va="center", fontsize=5.5, color=C.INK2)
ax.set_xticks(x)
ax.set_xticklabels([C.GENE_LABEL[g] for g in d.gene])
ax.set_ylim(0, 1)
ax.set_xlim(-0.4, len(d) - 0.3)
ax.set_ylabel("Global dN/dS (module-specific fit)")
from matplotlib.lines import Line2D  # noqa: E402
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C.INK2, label="channel"),
                   Line2D([], [], marker="o", ls="", color=C.BLUE, label="AQP4 extension"),
                   Line2D([], [], marker="o", ls="", color=C.ORANGE,
                          label="paralogue downstream frame")],
          loc="lower left", fontsize=5.5, handletextpad=0.1, bbox_to_anchor=(0, 1.0), ncol=2)
C.recessive_grid(ax)
C.save(fig, "Fig4", "c")
