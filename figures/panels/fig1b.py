"""Fig. 1b - Frequency of nested stop-codon contexts across set A2 (vertical bars)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
import numpy as np
import matplotlib.pyplot as plt

sp = C.species("aqp4")
a2 = sp[(sp.A2 == 1) & (sp.next4.str.len() == 4)]
n = len(a2)
ctx = a2["context"]
rows = [("UGA", ctx.str.startswith("UGA")),
        ("UGA-C", ctx.str.startswith("UGA-C")),
        ("UGA-CUA[G/U]", ctx.isin(["UGA-CUAG", "UGA-CUAU"])),
        ("UGA-CUAG", ctx.eq("UGA-CUAG")),
        ("UGA-CUAU", ctx.eq("UGA-CUAU"))]

fig, ax = plt.subplots(figsize=(2.4, 2.6))
x = np.arange(len(rows))
for i, (lab, m) in enumerate(rows):
    k = int(m.sum())
    pct = 100 * k / n
    col = C.ORANGE if lab == "UGA-CUAU" else C.BLUE
    ax.bar(i, pct, width=0.62, color=col)
    ax.text(i, pct + 2, f"{pct:.2f}%\n({k}/{n})", ha="center", va="bottom",
            fontsize=5.5, color=C.INK, linespacing=1.0)

ax.set_xticks(x)
ax.set_xticklabels([lab for lab, _ in rows], rotation=45, ha="right",
                   rotation_mode="anchor", fontsize=6.5, family="monospace")
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.6, len(rows) - 0.4)
ax.set_ylim(0, 125)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel(f"Species (% of A2, n = {n})")
C.recessive_grid(ax, "y")
C.save(fig, "Fig1", "b")
