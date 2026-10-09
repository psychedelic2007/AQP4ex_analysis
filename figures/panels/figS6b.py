"""Fig. S6b - Extension alignment (muscle) for representative species: human, one 28-codon-class
species per clade and two 11-codon acanthomorphs. All-gap columns among the shown rows are dropped."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

aln = C.alignment("aqp4", "muscle")
sp = C.species("aqp4")
sel = sp[sp.selection == 1]
reps = [C.REF]
for cl in C.CLADE_ORDER:
    d = sel[(sel.clade == cl) & (sel.species != C.REF)].sort_values("species")
    d28 = d[d.orf_codons == d.orf_codons.mode().iloc[0]] if len(d) else d
    if len(d28):
        reps.append(d28.species.iloc[0])
reps += [x for x in sel[sel.orf_codons == 11].sort_values("species").species if x not in reps][:2]
rows = [aln[s] for s in reps]
keep = [j for j in range(len(rows[0])) if any(r[j] != "-" for r in rows)]
fig, ax = plt.subplots(figsize=(0.13 * len(keep) + 1.8, 0.2 * len(reps) + 0.5))
for i, (s, r) in enumerate(zip(reps, rows)):
    for k, j in enumerate(keep):
        a = r[j]
        if a == "-":
            continue
        ax.add_patch(plt.Rectangle((k, -i - 0.45), 0.92, 0.9, color=C.AA_COLOR.get(a, C.NEUTRAL),
                                   alpha=0.25, lw=0))
        ax.text(k + 0.46, -i, a, ha="center", va="center", fontsize=5.5, family="monospace")
    lab = f"{s.replace('_', ' ')} ({C.CLADE_LABEL.get(sp.loc[s, 'clade'], '')}, {sp.loc[s, 'orf_codons']})"
    ax.text(-0.5, -i, lab, ha="right", va="center", fontsize=5.3,
            style="italic", fontweight="bold" if s == C.REF else "normal")
ax.set_xlim(-0.2, len(keep) + 0.2)
ax.set_ylim(-len(reps) + 0.4, 0.8)
ax.axis("off")
ax.set_title(f"MUSCLE: {len(aln[C.REF])} columns in full alignment, {len(keep)} shown",
             fontweight="normal", fontsize=6.5)
C.save(fig, "FigS6", "b")
