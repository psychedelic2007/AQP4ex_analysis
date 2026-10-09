"""Fig. S9b - Extension-only global dN/dS by phylogenetic depth (saturation check), with
profile-likelihood 95% CIs (results_final/aqp4/omega_ci/omega_profile_ci.tsv)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

pts = [("All (524)", "ext_mafft")] + [(f"{lab} ({C.fel_json('aqp4', f'depth_{k}_ext')['input']['number of sequences']})",
                                        f"depth_{k}_ext") for k, lab in C.DEPTH_LABEL.items()]
v = [C.global_omega("aqp4", f) for _, f in pts]
ci = pd.read_csv(os.path.join(C.RES, "aqp4", "omega_ci", "omega_profile_ci.tsv"), sep="\t",
                 index_col="fit")
lo = np.array([ci.loc[f, "ci95_lower"] for _, f in pts])
hi = np.array([ci.loc[f, "ci95_upper"] for _, f in pts])
assert np.allclose(v, [ci.loc[f, "omega_mle"] for _, f in pts], atol=1e-4), "CI table out of date"
fig, ax = plt.subplots(figsize=(2.8, 1.9))
ax.errorbar(range(len(v)), v, yerr=[np.array(v) - lo, hi - np.array(v)], color=C.BLUE,
            marker="o", ms=4, lw=1, elinewidth=0.8, capsize=2)
for i, y in enumerate(v):
    ax.text(i + 0.1, y, f"{y:.3f}", ha="left", va="center", fontsize=5.8)
ax.set_xticks(range(len(v)))
ax.set_xticklabels([p[0] for p in pts], rotation=30, ha="right", fontsize=5.8)
ax.set_xlim(-0.4, len(v) - 0.4)
ax.set_ylim(0.25, 0.5)
ax.set_ylabel("Extension dN/dS")
C.recessive_grid(ax)
C.save(fig, "FigS9", "b")
