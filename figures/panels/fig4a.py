"""Fig. 4a - Feature summary for AQP4 and five paralogues: UGA-C frequency (A2 with a
resolvable context) and frequency of the most common C-terminal hexapeptide (A4)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from collections import Counter  # noqa: E402
rows = []
for g in C.GENES:
    sp = C.species(g)
    a2 = sp[(sp.A2 == 1) & (sp.next4.str.len() == 4)]
    a4 = sp[sp.A4 == 1]
    ugac = int(a2.context.str.startswith("UGA-C").sum())
    top, k = Counter(s[-6:] for s in a4.extension if len(s) >= 6).most_common(1)[0]
    rows.append(dict(gene=g, ugac=100 * ugac / len(a2), n2=len(a2), k2=ugac,
                     hexa=100 * k / len(a4), top=top, n4=len(a4)))
d = pd.DataFrame(rows)
fig, axs = plt.subplots(1, 2, figsize=(4.4, 1.9), sharey=True)
for ax, col, title, lab in ((axs[0], "ugac", "UGA-C stop context", "% of species (A2)"),
                            (axs[1], "hexa", "Most frequent C-terminal hexapeptide",
                             "% of species (A4)")):
    y = np.arange(len(d))
    ax.barh(y, d[col], color=[C.BLUE if g == "aqp4" else C.NEUTRAL for g in d.gene],
            height=0.62)
    for i, r in d.iterrows():
        txt = (f"{r.ugac:.1f}%  ({r.k2}/{r.n2})" if col == "ugac"
               else f"{r.hexa:.1f}%  {r.top}  (n = {r.n4})")
        ax.text(r[col] + 1.5, i, txt, va="center", fontsize=5.5,
                family="monospace" if col == "hexa" else None)
    ax.set_xlim(0, 150)
    ax.set_xticks([0, 50, 100])
    ax.set_title(title, fontweight="normal", fontsize=6.5)
    ax.set_xlabel(lab)
    C.recessive_grid(ax, "x")
axs[0].set_yticks(range(len(d)))
axs[0].set_yticklabels([C.GENE_LABEL[g] for g in d.gene])
axs[0].invert_yaxis()
C.save(fig, "Fig4", "a")
