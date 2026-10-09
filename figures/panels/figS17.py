"""Fig. S17 - Share of tested downstream-frame columns under purifying selection, AQP4 against
the untranslated control frames of AQP1 and AQP2 (joint fits, occupancy >= 50%, q <= 0.20).
AQP4 is shown for the full selection set and at the depth of the paralogue sets (Tetrapoda,
Mammalia). Error bars are Wilson 95% intervals; brackets give two-sided Fisher's exact p."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

# (label, site tables pooled into the bar, is AQP4)
BARS = [("AQP4\nall", ["aqp4_joint_mafft_extension"], True),
        ("AQP4\nTetrapoda", ["aqp4_depth_tetrapods_joint_extension"], True),
        ("AQP4\nMammalia", ["aqp4_depth_mammals_joint_extension"], True),
        ("AQP1", ["aqp1_joint_mafft_extension"], False),
        ("AQP2", ["aqp2_joint_mafft_extension"], False),
        ("AQP1 +\nAQP2", ["aqp1_joint_mafft_extension", "aqp2_joint_mafft_extension"], False)]
# Comparisons quoted in the manuscript: (bar i, bar j)
TESTS = [(0, 5), (0, 3), (2, 5)]


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


rows = []
for lab, names, is4 in BARS:
    k = n = 0
    for nm in names:
        t = C.site_table(nm)
        k += int((t.tested & t.purifying).sum())
        n += int(t.tested.sum())
    lo, hi = wilson(k, n)
    rows.append(dict(lab=lab, k=k, n=n, p=100 * k / n, lo=100 * lo, hi=100 * hi, is4=is4))
    print(f"{lab.replace(chr(10), ' ')}: {k}/{n} = {100 * k / n:.1f}%")

fig, ax = plt.subplots(figsize=(3.6, 2.6))
x = np.arange(len(rows))
ax.bar(x, [r["p"] for r in rows], width=0.62, zorder=2,
       color=[C.BLUE if r["is4"] else C.NEUTRAL for r in rows])
ax.errorbar(x, [r["p"] for r in rows],
            yerr=[[r["p"] - r["lo"] for r in rows], [r["hi"] - r["p"] for r in rows]],
            fmt="none", ecolor=C.INK2, elinewidth=0.6, capsize=1.5, zorder=3)
for i, r in enumerate(rows):
    ax.text(i, r["hi"] + 1.5, f"{r['p']:.1f}%", ha="center", va="bottom", fontsize=6)
ax.set_xticks(x)
ax.set_xticklabels([f"{r['lab']}\n{r['k']}/{r['n']}" for r in rows])
ax.axvline(2.5, color=C.GRID, lw=0.6, zorder=1)

top = max(r["hi"] for r in rows) + 9
for level, (i, j) in enumerate(TESTS):
    a, b = rows[i], rows[j]
    p = fisher_exact([[a["k"], a["n"] - a["k"]], [b["k"], b["n"] - b["k"]]])[1]
    print(f"Fisher {a['lab'].replace(chr(10), ' ')} vs {b['lab'].replace(chr(10), ' ')}: p = {p:.3g}")
    y = top + level * 9
    ax.plot([i, i, j, j], [y - 1.5, y, y, y - 1.5], color=C.INK2, lw=0.6, zorder=3)
    ax.text((i + j) / 2, y + 0.8, f"p = {p:.2g}" if p >= 0.001 else f"p = {p:.1e}",
            ha="center", va="bottom", fontsize=6, color=C.INK)

ax.set_ylim(0, top + len(TESTS) * 9 + 2)
ax.set_yticks([0, 25, 50, 75, 100])
ax.spines["left"].set_bounds(0, 100)
ax.set_ylabel("Tested columns under\npurifying selection (%)")
ax.text(1, -0.25, "AQP4 extension", transform=ax.get_xaxis_transform(), ha="center",
        va="top", fontsize=6.5, color=C.INK2)
ax.text(4, -0.25, "Untranslated control frames", transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=6.5, color=C.INK2)
C.recessive_grid(ax)
C.save(fig, "FigS17", "")
