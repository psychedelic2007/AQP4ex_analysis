"""Fig. 2a - Length of the downstream reading frame (set A4), stacked by clade group,
with a kernel-density area curve (scaled to species counts) overlaid."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

sp = C.species("aqp4")
a4 = sp[sp.A4 == 1]
n = len(a4)
bins = np.arange(0, 46)

# Out-of-range values are silently dropped by the bins below - make that visible.
out = int((~a4.orf_codons.between(bins[0], bins[-1])).sum())
if out:
    print(f"[Fig2a] WARNING: {out}/{n} species have orf_codons outside "
          f"{bins[0]}-{bins[-1]} and are not shown")

fig, ax = plt.subplots(figsize=(3.5, 2.2))
bottom = np.zeros(len(bins))
for g in C.GROUP_ORDER:
    cnt = np.array([(a4[a4.group == g].orf_codons == b).sum() for b in bins])
    ax.bar(bins, cnt, bottom=bottom, width=0.8, color=C.GROUP_COLOR[g], label=g,
           edgecolor="white", linewidth=0.3, zorder=2)
    bottom += cnt

# --- area curve: Gaussian KDE in count units (bin width = 1 codon, so the curve
# is directly comparable to bar heights; its integral = number of species) ---
BW = 1.5  # bandwidth in codons; larger smooths more but flattens the peaks
x = a4.orf_codons.dropna().to_numpy(float)
grid = np.linspace(bins[0] - 0.5, bins[-1] + 0.5, 500)
dens = np.exp(-0.5 * ((grid[:, None] - x[None, :]) / BW) ** 2).sum(axis=1) \
       / (BW * np.sqrt(2 * np.pi))
ax.fill_between(grid, dens, color=C.INK, alpha=0.12, linewidth=0, zorder=1)
ax.plot(grid, dens, color=C.INK, linewidth=0.7, zorder=3)

for L, dy in ((28, 12), (11, 12), (43, 12)):
    k = int((a4.orf_codons == L).sum())
    ax.annotate(f"{L} codons\n{k} ({100 * k / n:.1f}%)", (L, bottom[L]), xytext=(0, dy),
                textcoords="offset points", ha="center", fontsize=5.8, color=C.INK)
ax.set_xlabel("Downstream reading-frame length (codons)")
ax.set_ylabel("Species")
ax.set_xlim(-1, 46)
ax.set_ylim(0, max(bottom.max(), dens.max()) * 1.2)
ax.tick_params(axis="both", labelsize=10)
ax.text(0.01, 0.97, f"n = {n}; median {np.median(a4.orf_codons):.0f}, mean "
        f"{a4.orf_codons.mean():.1f}", transform=ax.transAxes, fontsize=6, va="top")
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88), fontsize=5.5, handlelength=0.8)
C.recessive_grid(ax)
C.save(fig, "Fig2", "a")
