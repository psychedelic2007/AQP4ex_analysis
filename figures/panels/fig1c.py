"""Fig. 1c - Stop-codon context by clade (set A2), with n per clade."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import FancyBboxPatch, Rectangle

sp = C.species("aqp4")
a2 = sp[(sp.A2 == 1) & (sp.next4.str.len() == 4)]
cols = [("UGA-C", lambda c: c.str.startswith("UGA-C")),
        ("UGA-CUAG", lambda c: c.eq("UGA-CUAG")),
        ("UGA-CUAU", lambda c: c.eq("UGA-CUAU")),
        ("other", lambda c: ~c.isin(["UGA-CUAG", "UGA-CUAU"]))]
M, K, N = [], [], []
for cl in C.CLADE_ORDER:
    d = a2[a2.clade == cl]["context"]
    N.append(len(d))
    K.append([int(f(d).sum()) for _, f in cols])
    M.append([100 * k / len(d) if len(d) else np.nan for k in K[-1]])
M = np.array(M)

fig, ax = plt.subplots(figsize=(3.4, 3.1))
im = ax.imshow(M, cmap=C.SEQ_CMAP, vmin=0, vmax=100, aspect="auto")
ax.set_xticks(np.arange(-0.5, len(cols)), minor=True)
ax.set_yticks(np.arange(-0.5, len(N)), minor=True)
ax.grid(False, which="major")      
ax.grid(which="minor", color="black", linewidth=1.2)
ax.tick_params(which="minor", length=0)
for i in range(len(N)):
    for j in range(len(cols)):
        v = M[i, j]
        ax.text(j, i, f"{K[i][j]}/{N[i]}", ha="center", va="center", fontsize=5.5,
                color="white" if v > 60 else C.INK)
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([c for c, _ in cols], fontsize=6.5)
#ax.set_xticklabels([])
ax.xaxis.tick_top()
ax.set_yticks(range(len(N)))
ax.set_yticklabels([f"{C.CLADE_LABEL[c]} (n = {n})" for c, n in zip(C.CLADE_ORDER, N)])
#ax.set_yticklabels([])
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
cb.set_label("% of species in clade")
cb.outline.set_visible(True)
C.save(fig, "Fig1", "c")
