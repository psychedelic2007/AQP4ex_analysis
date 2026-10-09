"""Fig. S4b - Sequences of the 11-codon frames (all in Acanthomorphata)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from collections import Counter  # noqa: E402
sp = C.species("aqp4")
d = sp[(sp.A4 == 1) & (sp.orf_codons == 11)]
cnt = Counter(d.extension).most_common()
top = cnt[:10]
other = sum(c for _, c in cnt[10:])
labels = [k for k, _ in top] + ([f"{len(cnt) - 10} other sequences"] if other else [])
vals = [c for _, c in top] + ([other] if other else [])
fig, ax = plt.subplots(figsize=(3.2, 2.2))
ax.barh(range(len(vals)), vals, color=[C.MAGENTA] * len(top) + [C.NEUTRAL] * bool(other),
        height=0.65)
for i, v in enumerate(vals):
    ax.text(v + 0.8, i, str(v), va="center", fontsize=5.8)
ax.set_yticks(range(len(vals)))
ax.set_yticklabels(labels, family="monospace", fontsize=5.8)
ax.invert_yaxis()
ax.set_xlabel("Species")
clades = Counter(d.clade)
ax.set_title(f"11-codon frames: n = {len(d)}; " + ", ".join(f"{C.CLADE_LABEL[c]} {n}"
             for c, n in clades.items()), fontweight="normal", fontsize=6.2)
C.recessive_grid(ax, "x")
C.save(fig, "FigS4", "b")
