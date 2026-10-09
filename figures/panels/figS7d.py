"""Fig. S7d - global dN/dS of the AQP4 channel and extension from separate module-only
MG94xREV fits (same taxa and tree). Point estimates; the fits report no confidence interval.
Dashed line: global dN/dS of the joint (channel + extension) fit."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

vals = [("Channel\n(channel-only fit)", C.global_omega("aqp4", "channel"), C.INK2),
        ("Extension\n(extension-only fit)", C.global_omega("aqp4", "ext_mafft"), C.BLUE)]
joint = C.global_omega("aqp4", "joint_mafft")

fig, ax = plt.subplots(figsize=(2.2, 2.5))
C.recessive_grid(ax)
x = range(len(vals))
ax.bar(x, [v for _, v, _ in vals], width=0.55, color=[c for *_, c in vals], zorder=2)
for i, (_, v, _) in enumerate(vals):
    ax.text(i, v - 0.012, f"{v:.3f}", ha="center", va="top", fontsize=6, color="white")
ax.axhline(joint, ls="--", color=C.INK, lw=0.6, zorder=3)
ax.text(-0.47, joint + 0.008, f"joint fit {joint:.3f}", ha="left", va="bottom", fontsize=5.5)
ax.set_xticks(list(x))
ax.set_xticklabels([l for l, *_ in vals], fontsize=6)
ax.set_xlim(-0.5, 1.5)
ax.set_ylim(0, 0.5)
ax.set_ylabel("Global dN/dS (ω)")
ax.set_title(f"AQP4 module dN/dS (extension/channel = {vals[1][1] / vals[0][1]:.1f}×)",
             fontweight="normal", fontsize=6)
C.save(fig, "FigS7", "d")
