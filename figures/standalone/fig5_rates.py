#!/usr/bin/env python3
"""Fig. 5a-f (left panels) - site-wise synonymous (alpha, open circles) and non-synonymous
(beta, filled orange circles) rates at the tested downstream-frame codons of each gene,
in human residue numbering. The structure panels on the right are the AlphaFold 3 models
coloured in PyMOL with figures/pymol/<gene>_{channel,tail}_{alpha,beta}.pml.

Input: figures/pymol/<gene>_tail_alpha_sites.tsv, written by scripts/fel_to_pymol.py from
the joint FEL fit of record. Only tested codons are drawn; columns that are gaps in the
human sequence have no human residue number and are therefore not shown here.

    python3 figures/standalone/fig5_rates.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PYMOL = os.path.join(ROOT, "figures", "pymol")
OUT = os.path.join(ROOT, "figures", "output", "standalone")
PAPER = os.path.join(ROOT, "figures", "paper")
GENES = ["aqp1", "aqp2", "aqp4", "aqp5", "aqp6", "aqp8"]          # rows a-f of Fig. 5

os.makedirs(OUT, exist_ok=True)
os.makedirs(PAPER, exist_ok=True)
for letter, g in zip("abcdef", GENES):
    path = os.path.join(PYMOL, f"{g}_tail_alpha_sites.tsv")
    if not os.path.exists(path):
        sys.exit(f"missing {path}; run scripts/fel_to_pymol.py first")
    df = pd.read_csv(path, sep="\t")
    df = df[df["call"] != "not tested"].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.5, 3.0))
    for i, r in df.iterrows():
        ax.plot([i, i], [r.alpha, r.beta], color="#E0E0E0", lw=1.5, zorder=1)
    ax.scatter(df.index, df.alpha, facecolors="none", edgecolors="#404040", lw=0.9, s=28,
               label="α (synonymous)", zorder=2)
    ax.scatter(df.index, df.beta, color="#E6613B", s=28, label="β (non-synonymous)", zorder=2)
    ax.axhline(1.0, color="black", ls="--", lw=0.6, zorder=0)
    ax.set_yscale("log")
    ax.set_xlim(-1, len(df))
    step = max(1, len(df) // 20)
    ax.set_xticks(range(0, len(df), step))
    ax.set_xticklabels(df.human_residue.iloc[::step], rotation=90, fontsize=7)
    ax.set_ylabel("Rate")
    ax.set_title(f"{g.upper()}   α = {df.alpha.min():.2f}–{df.alpha.max():.2f}   "
                 f"β = {df.beta.min():.2f}–{df.beta.max():.2f}", loc="left", fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("svg", "pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{g}_tail_alpha_beta.{ext}"), dpi=300)
        fig.savefig(os.path.join(PAPER, f"Fig5{letter}_rates.{ext}"), dpi=300)
    plt.close(fig)
    print(f"saved figures/paper/Fig5{letter}_rates.{{svg,pdf,png}}  ({g}, {len(df)} tested codons)")
