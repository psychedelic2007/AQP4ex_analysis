"""Fig. 2c - Sequence of the extension, anchored at the C-terminus (no alignment needed):
probability logos for the 28-codon and 11-codon classes, and the most frequent
C-terminal hexapeptides and N-terminal pentapeptides (denominator: A4)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import logomaker  # noqa: E402
sp = C.species("aqp4")
a4 = sp[sp.A4 == 1]
n = len(a4)
AA = "ACDEFGHIKLMNPQRSTVWY"
NSPLIT = 5

def prob(seqs, L):
    m = np.zeros((L, len(AA)))
    for s in seqs:
        for i, a in enumerate(s[-L:]):
            if a in AA:
                m[i, AA.index(a)] += 1
    return pd.DataFrame(m / m.sum(1, keepdims=True), columns=list(AA))


fig = plt.figure(figsize=(8.1, 2.7))
gs = fig.add_gridspec(2, 3, width_ratios=[2.8, 1.2, 1.0], height_ratios=[1, 1], wspace=0.35,
                      hspace=0.75)
for row, L in enumerate((28, 11)):
    seqs = a4[a4.orf_codons == L].extension
    ax = fig.add_subplot(gs[row, 0] if L == 28 else gs[row, 0])
    df = prob(seqs, L)
    if L == 11:   # right-align the 11-codon logo under the last 11 positions of the 28-codon logo
        full = pd.DataFrame(0.0, index=np.arange(28), columns=df.columns)
        full.iloc[:NSPLIT] = df.iloc[:NSPLIT].values
        full.iloc[28 - (L - NSPLIT):] = df.iloc[NSPLIT:].values
        df = full
    logomaker.Logo(df, ax=ax, color_scheme=C.AA_COLOR, vpad=0.05, width=0.85)
    if L == 11:
        mid = (NSPLIT - 1 + 28 - (L - NSPLIT)) / 2
        ax.text(mid, 0.5, "· · ·", ha="center", va="center", fontsize=8, color=C.NEUTRAL)
    ax.set_xlim(-0.6, 27.6)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Probability")
    ax.set_yticks([0, 0.5, 1])
    if L == 28:
        ax.set_xticks(range(0, 28, 3))
        ax.set_xticklabels([325 + i for i in range(0, 28, 3)], fontsize=10)
        #ax.set_xticklabels([])
        ax.set_xlabel("Human position", labelpad=1)
    else:
        cpos = list(range(28 - (L - NSPLIT), 28))
        ax.set_xticks(list(range(NSPLIT)) + cpos)
        ax.set_xticklabels(["N"] + [f"+{i}" for i in range(1, NSPLIT)]
                           + [f"−{27 - i}" if i < 27 else "C" for i in cpos], fontsize=10, rotation=45)
        #ax.set_xticklabels([])
        ax.set_xlabel("Residues from C-terminus", labelpad=1)
    ax.set_title(f"{L}-codon frames (n = {len(seqs)}, {100 * len(seqs) / n:.1f}% of A4)",
                 fontweight="normal", fontsize=6.5)


def topbars(ax, counts, title, hl):
    top = counts.most_common(6)[::-1]
    ax.barh(range(len(top)), [100 * c / n for _, c in top],
            color=[C.BLUE if k in hl else C.NEUTRAL for k, _ in top], height=0.65)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([k for k, _ in top], family="monospace", fontsize=6)
    #ax.set_yticklabels([])
    for i, (_, c) in enumerate(top):
        ax.text(100 * c / n + 1, i, f"{100 * c / n:.1f}%", va="center", fontsize=5.5)
    ax.set_title(title, fontweight="normal", fontsize=6.5)
    ax.set_xlabel(f"% of A4 (n = {n})")
    ax.set_xlim(0, 95)
    ax.tick_params(axis="both", labelsize=10)


from collections import Counter  # noqa: E402
hexa = Counter(s[-6:] for s in a4.extension if len(s) >= 6)
penta = Counter(s[:5] for s in a4.extension if len(s) >= 5)
topbars(fig.add_subplot(gs[:, 1]), hexa, "C-terminal hexapeptide", {"KETDLL"})
topbars(fig.add_subplot(gs[:, 2]), penta, "N-terminal pentapeptide", {"LEDRT"})
C.save(fig, "Fig2", "c")
