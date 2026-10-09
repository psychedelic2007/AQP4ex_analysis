"""Fig. S15b - Probability logos of the stop codon and the following nucleotide (+4) in
species whose stop codon is UGA, per gene (set A2, resolvable stop context). Letter
height = fraction of species, so the +4 letters read directly as the percentages
printed above them. Positions +1..+3 are 1.0 by construction (only UGA species).
Set MODE = "bits" for the information-content version (letter height = p x R, with the
Schneider small-sample correction e(n) = 3 / (2 ln2 n)); there a mixed +4 position is
drawn short even when its commonest base is frequent."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import logomaker  # noqa: E402

ORDER = ["aqp1", "aqp2", "aqp4", "aqp5", "aqp6", "aqp8"]
MODE = "probability"   # or "bits"
NT = "ACGU"
NT_COLOR = {"A": C.GREEN, "C": C.BLUE, "G": C.YELLOW, "U": C.RED}


def info(seqs):
    n = len(seqs)
    cnt = np.array([[sum(s[i] == b for s in seqs) for b in NT] for i in range(4)], float)
    p = cnt / n
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -np.nansum(np.where(p > 0, p * np.log2(p), 0), axis=1)
    R = np.clip(2 - (H + 3 / (2 * np.log(2) * n)), 0, None)
    h = p if MODE == "probability" else p * R[:, None]
    return pd.DataFrame(h, columns=list(NT)), p


TOP = 1.0 if MODE == "probability" else 2.0
fig, axes = plt.subplots(1, 6, figsize=(7.1, 1.7), sharey=True)
for ax, g in zip(axes, ORDER):
    sp = C.species(g)
    ctx = sp[(sp.A2 == 1) & (sp.context != "")].context
    seqs = [c[:3] + c[4] for c in ctx if c.startswith("UGA")]
    df, p = info(seqs)
    logomaker.Logo(df, ax=ax, color_scheme=NT_COLOR, vpad=0.05, width=0.85)
    top = int(p[3].argmax())
    print(f"{g.upper()}: n={len(seqs)} +4 " +
          " ".join(f"{b}={100 * p[3, i]:.1f}%" for i, b in enumerate(NT)))
    ax.set_title(f"{C.GENE_LABEL[g]}  (n = {len(seqs)})", fontsize=6.5)
    ax.text(3, TOP * 1.06, f"{NT[top]} {100 * p[3, top]:.1f}%", ha="center", fontsize=5.5,
            color=C.INK2)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["+1", "+2", "+3", "+4"])
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(0, TOP * 1.175)
    ax.set_yticks([0, 0.5, 1] if MODE == "probability" else [0, 1, 2])
axes[0].set_ylabel("Probability" if MODE == "probability" else "Information (bits)")
fig.supxlabel("Position (stop codon +1 to +3, next nucleotide +4)", fontsize=7, y=0.02)
fig.tight_layout(w_pad=0.6)
C.save(fig, "FigS15", "b")
