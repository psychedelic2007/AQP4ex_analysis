"""Fig. 3b - Residue composition at the 28 tested columns (MAFFT alignment, 524 species,
gaps excluded); human residue shown below, open circle where it is not the dominant state."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import logomaker  # noqa: E402
t = C.site_table("aqp4_joint_mafft_extension")
t = t[t.tested]
cols = (t.column.astype(int) - 1).tolist()
aln = C.alignment("aqp4", "mafft")
AA = "ACDEFGHIKLMNPQRSTVWY"
M = np.zeros((len(cols), len(AA)))
for s in aln.values():
    for i, c in enumerate(cols):
        if s[c] in AA:
            M[i, AA.index(s[c])] += 1
occ = M.sum(1)
df = pd.DataFrame(M / occ[:, None], columns=list(AA))

fig, ax = plt.subplots(figsize=(7.1, 1.9))
logomaker.Logo(df, ax=ax, color_scheme=C.AA_COLOR, vpad=0.05, width=0.85)
ax.set_ylim(0, 1)
ax.set_ylabel("Frequency")
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([f"{a}\n{p}" for a, p in zip(t.human_aa, t.human_residue)], fontsize=5.8)
ax.set_xlim(-0.6, len(cols) - 0.4)
for i, (_, r) in enumerate(t.iterrows()):
    if r.dominant_aa != r.human_aa:
        ax.scatter([i], [-0.30], s=14, facecolor="white", edgecolor=C.ORANGE, lw=0.9,
                   clip_on=False, transform=ax.get_xaxis_transform())
ax.text(-1.0, -0.30, "human ≠\ndominant", transform=ax.get_xaxis_transform(), ha="right",
        va="center", fontsize=5.3, color=C.INK2)
ax.set_title(f"n = {len(aln)} species; occupancy of tested columns {int(occ.min())}–{int(occ.max())}",
             fontweight="normal", fontsize=6.3)
C.save(fig, "Fig3", "b")
