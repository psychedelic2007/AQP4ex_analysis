"""Fig. S11b - Downstream reading-frame length distributions (A4) for each gene."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

fig, axs = plt.subplots(2, 3, figsize=(5.4, 3.0), sharex=True)
for ax, g in zip(axs.flat, C.GENES):
    a4 = C.species(g).query("A4 == 1")
    ax.hist(a4.orf_codons.clip(upper=100), bins=np.arange(-0.5, 101.5, 1),
            color=C.BLUE if g == "aqp4" else C.ORANGE)
    ax.set_title(f"{C.GENE_LABEL[g]} (n = {len(a4)}, median {a4.orf_codons.median():.0f})",
                 fontweight="normal", fontsize=6.3)
    ax.tick_params(labelsize=5.5)
for ax in axs[1]:
    ax.set_xlabel("Codons (>100 pooled)")
fig.supylabel("Species", fontsize=6.5)
fig.tight_layout(h_pad=0.6)
C.save(fig, "FigS11", "b")
