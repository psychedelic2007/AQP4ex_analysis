"""Fig. S7a - Positive control: site-wise alpha and beta along the AQP4 channel domain
(joint fit), with calls and the two NPA motifs of the human sequence."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

t = C.site_table("aqp4_joint_mafft_channel")
t["pos"] = t.human_residue.astype(int)
human = "".join(t.human_aa)
npa = [i + 1 for i in range(len(human) - 2) if human[i:i + 3] == "NPA"]
fig, ax = plt.subplots(figsize=(7.1, 2.1))
tt = t[t.tested]
ax.scatter(tt.pos, tt.alpha, s=4, facecolor="white", edgecolor=C.INK2, lw=0.5, label="α")
ax.scatter(tt.pos, tt.beta, s=4, color=C.ORANGE, lw=0, label="β")
for _, r in tt.iterrows():
    if r.purifying:
        ax.plot([r.pos, r.pos], [25, 40], color=C.BLUE if r.bonf else "#9ec5f4", lw=0.8)
    elif r.diversifying:
        ax.plot([r.pos, r.pos], [25, 40], color=C.RED, lw=1.4)
        ax.annotate(f"codon {r.pos} diversifying (q = {r.q:.3f})", (r.pos, 40), xytext=(r.pos + 5, 60),
                    fontsize=5.5, color=C.RED)
for p in npa:
    ax.axvspan(p - 0.5, p + 2.5, color=C.YELLOW, alpha=0.3, lw=0)
    ax.text(p + 1, 0.0004, "NPA", ha="center", fontsize=5.5, color=C.INK2)
ax.set_yscale("log")
ax.set_ylim(0.0002, 90)
ax.set_xlim(0, 324)
ax.set_xlabel("Human AQP4 residue")
ax.set_ylabel("Rate")
n_t, n_p = int(t.tested.sum()), int(t.purifying.sum())
ax.set_title(f"Tested {n_t}; purifying {n_p} ({100 * n_p / n_t:.1f}%), Bonferroni "
             f"{int((t.purifying & t.bonf).sum())}; untested columns omitted. Bars: dark = "
             "Bonferroni, light = q ≤ 0.20", fontweight="normal", fontsize=6)
ax.legend(loc="lower right", ncol=2, fontsize=5.5, markerscale=2)
C.save(fig, "FigS7", "a")
