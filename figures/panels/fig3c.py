"""Fig. 3c - Synonymous rates of tested channel vs extension codons in the same joint fit."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from scipy.stats import mannwhitneyu  # noqa: E402
ch = C.site_table("aqp4_joint_mafft_channel")
ex = C.site_table("aqp4_joint_mafft_extension")
a_ch = ch[ch.tested].alpha.values
a_ex = ex[ex.tested].alpha.values
p_mw = mannwhitneyu(a_ex, a_ch, alternative="two-sided").pvalue
rng = np.random.default_rng(20260926)
obs = np.median(a_ch) - np.median(a_ex)
pool = np.concatenate([a_ch, a_ex])
k = 0
for _ in range(10000):
    rng.shuffle(pool)
    if abs(np.median(pool[:len(a_ch)]) - np.median(pool[len(a_ch):])) >= abs(obs):
        k += 1
p_perm = (k + 1) / 10001

fig, ax = plt.subplots(figsize=(2.3, 2.4))
for i, (v, col, lab) in enumerate(((a_ch, C.INK2, "channel"), (a_ex, C.BLUE, "extension"))):
    vp = ax.violinplot([v], positions=[i], widths=0.75, showextrema=False)
    for b in vp["bodies"]:
        b.set_facecolor(col)
        b.set_alpha(0.15)
        b.set_edgecolor("none")
    ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=5 if i == 0 else 10, color=col,
               alpha=0.5 if i == 0 else 0.85, lw=0)
    ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=C.INK, lw=1.4)
    ax.text(i + 0.3, np.median(v), f"{np.median(v):.3f}", va="center", fontsize=5.8)
r325 = ex[ex.human_residue.astype(str) == "325"].iloc[0]
ax.annotate("L325 (α = 0.029)", (1, r325.alpha), xytext=(1.0, 0.012), fontsize=5.5,
            ha="center", color=C.INK2)
ax.set_yscale("log")
ax.set_xticks([0, 1])
ax.set_xticklabels([f"channel\n(n = {len(a_ch)})", f"extension\n(n = {len(a_ex)})"])
ax.set_ylabel("Synonymous rate α")
ax.set_xlim(-0.6, 1.8)
ax.set_title(f"Mann–Whitney P = {p_mw:.1e}\npermutation P = {p_perm:.4f}",
             fontweight="normal", fontsize=6.2)
C.recessive_grid(ax)
C.save(fig, "Fig3", "c")
