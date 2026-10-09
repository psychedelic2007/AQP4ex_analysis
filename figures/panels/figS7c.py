"""Fig. S7c - alpha vs beta for all tested codons of the AQP4 joint (channel + extension) fit.
Channel codons as circles, extension codons as diamonds, coloured by call. The global dN/dS
in the title is the joint model's MG94xREV parameter, not a summary of the plotted points."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

FLOOR = 1e-4
mods = {m: C.site_table(f"aqp4_joint_mafft_{m}") for m in ("channel", "extension")}
tested = {m: t[t.tested] for m, t in mods.items()}
calls = (("not called", C.NEUTRAL, lambda t: ~t.purifying & ~t.diversifying),
         ("purifying", C.BLUE, lambda t: t.purifying),
         ("diversifying", C.RED, lambda t: t.diversifying))

fig, ax = plt.subplots(figsize=(3.0, 2.8))
for _, col, sel in calls:
    t = tested["channel"]
    ax.scatter(t[sel(t)].alpha, t[sel(t)].beta.clip(lower=FLOOR), s=5, color=col, lw=0,
               alpha=0.6, zorder=2)
for _, col, sel in calls:
    t = tested["extension"]
    ax.scatter(t[sel(t)].alpha, t[sel(t)].beta.clip(lower=FLOOR), s=16, marker="D",
               facecolor=col, edgecolor=C.INK, lw=0.5, zorder=3)
ax.plot([FLOOR, 20], [FLOOR, 20], ls="--", color=C.INK2, lw=0.5, zorder=1)
ax.text(1.5e-3, 2.2e-3, "ω = 1", rotation=45, rotation_mode="anchor", fontsize=5.5, color=C.INK2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("α (synonymous)")
ax.set_ylabel(r"β (non-synonymous; floored at $\mathregular{10^{-4}}$)")

handles = []
for m, marker, size in (("channel", "o", 3.5), ("extension", "D", 4)):
    t = tested[m]
    n_tot, n_p, n_d = len(mods[m]), int(t.purifying.sum()), int(t.diversifying.sum())
    handles.append(Line2D([], [], ls="", marker=marker, ms=size, mfc="white",
                          mec=C.INK if m == "extension" else C.INK2, mew=0.6,
                          label=f"{m} ({len(t)}/{n_tot} tested)\n{n_p} purifying, {n_d} diversifying"))
for lab, col, _ in calls:
    handles.append(Line2D([], [], ls="", marker="s", ms=4, color=col, label=lab))
ax.legend(handles=handles, fontsize=5, loc="upper left", handletextpad=0.3, frameon=False,
          borderaxespad=0.3, labelspacing=0.6)

n_floor = sum(int((t.beta < FLOOR).sum()) for t in tested.values())
ax.set_title(f"AQP4 joint fit, {sum(len(t) for t in tested.values())} tested codons\n"
             f"Model global dN/dS = {C.global_omega('aqp4', 'joint_mafft'):.3f} (fitted parameter, "
             f"not computed from points)\n{n_floor} codons with "
             r"β < $\mathregular{10^{-4}}$ plotted at the floor",
             fontweight="normal", fontsize=5.5)
C.save(fig, "FigS7", "c")
