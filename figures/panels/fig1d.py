"""Fig. 1d - UGA-CUAG on the time-calibrated species tree: the eight losses common to
all most-parsimonious reconstructions (filled) and the unresolved ninth change
(hollow: gain on the Gnathostomata stem OR loss on the lamprey stem)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

t, state = C.cuag_tree()
cost = C.sankoff(t, state)
ev_present = C.parsimony_events(t, state, 1, cost)     # 9 losses
ev_absent = C.parsimony_events(t, state, 0, cost)      # 1 gain + 8 losses
key = lambda ev: {(frozenset(x.name for x in c.get_terminals()), k) for c, k in ev}
shared = key(ev_present) & key(ev_absent)
common = [(c, k) for c, k in ev_present
          if (frozenset(x.name for x in c.get_terminals()), k) in shared]
alt = [(c, k) for c, k in ev_present + ev_absent
       if (frozenset(x.name for x in c.get_terminals()), k) not in shared]

fig, ax = plt.subplots(figsize=(3.6, 5.6))
OFFSET = {"Latimeria_chalumnae": ((0, 7), "center"),
          "Petromyzontidae": ((-6, 0), "right"),
          "Rattus_rattus": ((-6, 0), "right"),
          "Rhinoderma_darwinii": ((-6, 4), "right"),
          "Alligator_mississippiensis": ((-6, 0), "right"),
          "Gnathostomata": ((-6, 0), "right")}
pos, x, y = C.draw_tree(ax, t, state, events=common, ambiguous=alt)
C.clade_strip(ax, t, max(x.values()) + 12, width=6)
ax.set_xlim(-5, max(x.values()) + 18)
for ch, kind, xm, ym in pos:
    tips = [tp.name for tp in ch.get_terminals()]
    lab = tips[0].replace("_", " ") if len(tips) == 1 else f"{C.lca_name(tips)} ({len(tips)})"
    off, ha = OFFSET.get(tips[0] if len(tips) == 1 else C.lca_name(tips), ((-5, -5), "right"))
    ax.annotate(lab, (xm, ym), xytext=off, textcoords="offset points", ha=ha, va="center",
                fontsize=5.2, color=C.INK2, style="italic" if len(tips) == 1 else "normal",
                bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.85))
from matplotlib.lines import Line2D  # noqa: E402
h = [Line2D([], [], marker="s", ls="", color=C.BLUE, ms=3, label="UGA-CUAG present"),
     Line2D([], [], marker="s", ls="", color=C.ORANGE, ms=3, label="UGA-CUAG absent"),
     Line2D([], [], marker="v", ls="", mfc=C.ORANGE, mec="white", ms=5,
            label="loss in all MP reconstructions"),
     Line2D([], [], marker="^", ls="", mfc="white", mec=C.BLUE, ms=5,
            label="gain (root absent) or"),
     Line2D([], [], marker="v", ls="", mfc="white", mec=C.ORANGE, ms=5,
            label="loss (root present): equal cost")]
ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.5,
          handletextpad=0.2, columnspacing=0.8)
ax.text(0, 1.105, f"{len(state)} tips; parsimony cost 9 under both root states",
        transform=ax.transAxes, fontsize=5.5, color=C.INK2)
C.save(fig, "Fig1", "d")
