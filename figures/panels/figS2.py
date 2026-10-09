"""Fig. S2 - The five AQP4 records lacking UGA-C, aligned at the stop codon with the most
closely related A2 species that carries UGA-C (curated records preferred)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

sp = C.species("aqp4")
a2 = sp[(sp.A2 == 1) & (sp.next4.str.len() == 4)]
bad = a2[~a2.context.str.startswith("UGA-C")].species.tolist()
good = a2[a2.context.str.startswith("UGA-C")]
cds, utr = C.raw_cds("aqp4"), C.raw_utr("aqp4")
lin = {s: C.lineage(s) for s in set(bad) | set(good.species)}


def shared(a, b):
    n = 0
    for x, y in zip(lin[a], lin[b]):
        if x != y:
            break
        n += 1
    return n


pairs = []
for s in bad:
    best = max(good.species, key=lambda o: (shared(s, o), good.loc[o, "accession"].startswith("NM_")))
    pairs.append((s, best))


def window(s):
    acc = sp.loc[s, "accession"]
    return cds[acc][-15:-3], cds[acc][-3:], utr[acc][:12], acc


fig, ax = plt.subplots(figsize=(6.2, 0.42 * 2 * len(pairs) + 0.4))
ax.axis("off")
y = 0
for s, o in pairs:
    for k, name in enumerate((s, o)):
        up, stop, down, acc = window(name)
        yy = -y
        ax.text(0.0, yy, name.replace("_", " "), fontsize=6.2, style="italic", va="center",
                color=C.INK if k == 0 else C.INK2)
        ax.text(0.30, yy, acc, fontsize=5.8, va="center", color=C.INK2)
        seq = up + stop + down
        for i, nt in enumerate(seq):
            x = 0.46 + i * 0.0175
            in_stop = 12 <= i < 15
            ctx = 15 <= i < 19
            ref = window(o)[0] + window(o)[1] + window(o)[2]
            diff = k == 0 and nt != ref[i]
            if in_stop or ctx:
                ax.add_patch(plt.Rectangle((x - 0.008, yy - 0.3), 0.0165, 0.6,
                                           color=C.ORANGE if in_stop else "#cde2fb", lw=0,
                                           alpha=0.35 if in_stop else 1))
            ax.text(x, yy, nt.replace("T", "U"), family="monospace", fontsize=6.3, ha="center",
                    va="center", color=C.RED if diff else C.INK,
                    fontweight="bold" if diff else "normal")
        y += 0.8
    y += 0.5
ax.text(0.46 + 13 * 0.0175, 0.8, "stop", ha="center", fontsize=5.8, color=C.INK2)
ax.text(0.46 + 16.5 * 0.0175, 0.8, "+4 to +7", ha="center", fontsize=5.8, color=C.INK2)
ax.set_ylim(-y + 0.2, 1.2)
ax.set_xlim(0, 1)
ax.text(0, -y - 0.1, "Each record lacking UGA-C (black) is followed by its closest relative with UGA-C "
        "(grey); differences in red.", fontsize=5.5, color=C.INK2)
C.save(fig, "FigS2", "")
