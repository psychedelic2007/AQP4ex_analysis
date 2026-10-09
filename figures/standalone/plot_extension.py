#!/usr/bin/env python3
"""
Figure plots for Result 2 — extension length and terminal sequence conservation.

Companion to plot_stopcontext.py, same conventions: ONE FILE PER PLOT, no
panels, sized for assembly in Inkscape, plus a text file recording every
plotted number.

Usage
-----
    python3 plot_extension.py --gene aqp4
    python3 plot_extension.py --gene aqp4 --datadir /path/to/fasta
    python3 plot_extension.py --gene aqp4 --only length,logo_cterm
    python3 plot_extension.py --gene aqp4 --compare aqp1 aqp2 aqp5 aqp6 aqp8

Plots produced
--------------
  MAIN
    length          length distribution against the composition null
    length_ranked   ranked frequency of the commonest discrete lengths
    logo_cterm      sequence logo, terminal 8 residues
    logo_nterm      sequence logo, first 8 residues
    asymmetry       commonest C- vs N-terminal motif frequency
    motif_cterm     ranked C-terminal hexapeptide frequency

  SUPPLEMENTARY
    null_vs_gc      chance expectation against 3'UTR GC, all genes
    length_ecdf     cumulative length distribution
    motif_nterm     ranked N-terminal pentapeptide frequency
    conservation    per-position conservation along the extension

Denominator
-----------
All plots use A4: species with a CDS, a 3'UTR record, a reading frame ending
at a genuine in-frame stop codon, and no frameshift. Stated in every legend.
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit("biopython required:  pip install biopython")

ACC = re.compile(r"([A-Za-z]{2}_[\d.]+)")
ORG = re.compile(r"\[organism=([^\]]+)\]")
REF_LEN = dict(aqp1=269, aqp2=271, aqp4=323, aqp5=265, aqp6=282, aqp8=261)
LONG_CUT = 35          # reading frames above this are inspected individually

# ---------------------------------------------------------------- styling --
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 150,
})

C = dict(hi="#1b7f79", lo="#9aa5b1", alt="#8b3a62", warn="#c05746",
         null="#5d708f", nterm="#c77b3a", cterm="#1b7f79")

# Amino acid colouring by physicochemical class, for the logos.
AA_COL = {}
for grp, col in (("AVLIMFWP", "#1b7f79"),      # hydrophobic
                 ("STNQCGY", "#5d708f"),        # polar
                 ("DE", "#c53030"),             # acidic
                 ("KRH", "#2b6cb0")):           # basic
    for a in grp:
        AA_COL[a] = col
AA20 = "ACDEFGHIKLMNPQRSTVWY"


def sanitise(s):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", s)).strip("_")


def save(fig, outdir, gene, name):
    for ext in ("svg", "pdf"):
        fig.savefig(os.path.join(outdir, f"{gene}_{name}.{ext}"),
                    bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  [ok] {gene}_{name}.svg / .pdf")


# ------------------------------------------------------------------ loader --
def load(gene, datadir):
    """species -> extension sequence, length, 3'UTR GC, QC flags.

    Returns the full set; filtering to A4 happens in a4().
    """
    utr = {}
    for r in SeqIO.parse(os.path.join(datadir, f"3p_utr_{gene}.fna"), "fasta"):
        a = ACC.search(r.id) or ACC.search(r.description)
        if a:
            utr[a.group(1)] = str(r.seq).upper()

    valid = {}
    for r in SeqIO.parse(os.path.join(datadir, f"cds_{gene}.fna"), "fasta"):
        a, o = ACC.search(r.id), ORG.search(r.description)
        if not a or not o:
            continue
        s = str(r.seq).upper()
        if len(s) % 3 or not s.startswith("ATG"):
            continue
        p = str(Seq(s).translate())
        if not p.endswith("*") or "*" in p[:-1]:
            continue
        valid[a.group(1)] = (sanitise(o.group(1)), s, p[:-1])

    reflen = REF_LEN.get(gene, 323)
    best = {}
    for acc, (sp, s, p) in valid.items():
        k = (not acc.startswith("NM_"), abs(len(p) - reflen))
        if sp not in best or k < best[sp][0]:
            best[sp] = (k, acc, s, p)

    out = {}
    for sp, (_, acc, s, p) in best.items():
        u = utr.get(acc, "")
        if not u:
            continue
        tail = u[:len(u) // 3 * 3]
        ext = str(Seq(tail).translate(to_stop=True)) if tail else ""
        used = 3 * len(ext)
        stop_ok = (len(tail) >= used + 3
                   and str(Seq(tail[used:used + 3]).translate()) == "*")
        gc = ((u[:100].count("G") + u[:100].count("C")) / len(u[:100])
              if len(u) >= 100 else None)
        out[sp] = dict(acc=acc, ext=ext, n=len(ext), stop_ok=stop_ok,
                       gc=gc, utr_len=len(u))
    return out


def a4(rows):
    """Apply the A4 filter: confirmed stop, no frameshift."""
    motif = re.compile(r"[A-Z]ET[DEN][LI][LI]?$")
    keep, dropped = {}, []
    for sp, r in rows.items():
        if not r["stop_ok"]:
            continue
        if r["n"] > LONG_CUT and not motif.search(r["ext"]):
            dropped.append((sp, r["n"]))
            continue
        keep[sp] = r
    return keep, dropped


def comp_null(rows):
    """Expected codons before an in-frame stop, from 3'UTR base composition."""
    gc = [r["gc"] for r in rows.values() if r["gc"] is not None]
    if not gc:
        return None, None, None
    g = float(np.mean(gc))
    a = t = (1 - g) / 2
    c = gg = g / 2
    p = t * a * a + t * a * gg + t * gg * a
    return 1 / p, g, p


# -------------------------------------------------------------------- logo --
def letter(ch, x, y, w, h, ax, colour):
    tp = TextPath((0, 0), ch, size=1,
                  prop=matplotlib.font_manager.FontProperties(
                      family="DejaVu Sans", weight="bold"))
    bb = tp.get_extents()
    if bb.width == 0 or bb.height == 0:
        return
    t = (Affine2D().translate(-bb.x0, -bb.y0)
         .scale(w / bb.width, h / bb.height).translate(x, y))
    ax.add_patch(PathPatch(t.transform_path(tp), facecolor=colour,
                           edgecolor="none"))


def aa_logo(seqs, labels, outdir, gene, name, rep, note, width=4.2):
    """Information-content logo over aligned amino acid columns.

    Maximum information for 20 amino acids is log2(20) = 4.32 bits, but real
    protein columns rarely approach it; the axis is scaled to the data.
    """
    n, L = len(seqs), len(seqs[0])
    fig, ax = plt.subplots(figsize=(width, 2.3))
    lines, peak = [], 0
    for i in range(L):
        col = Counter(s[i] for s in seqs if s[i] in AA20)
        tot = sum(col.values())
        if tot == 0:
            lines.append(f"pos {labels[i]:<8} no data")
            continue
        freqs = {a: col.get(a, 0) / tot for a in AA20}
        H = -sum(p * np.log2(p) for p in freqs.values() if p > 0)
        R = max(np.log2(20) - H - (19 / (2 * np.log(2) * tot)), 0)
        peak = max(peak, R)
        y = 0
        for a, p in sorted(freqs.items(), key=lambda t: t[1]):
            if p <= 0.01:
                continue
            h = p * R
            letter(a, i + 0.12, y, 0.76, h, ax, AA_COL.get(a, "#666"))
            y += h
        top = col.most_common(3)
        lines.append(f"pos {labels[i]:<8} bits={R:5.3f}  n={tot:<5} " +
                     "  ".join(f"{a}={100*c/tot:.1f}%" for a, c in top))

    ax.set_xlim(0, L)
    ax.set_ylim(0, max(peak * 1.1, 0.5))
    ax.set_xticks(np.arange(L) + 0.5)
    ax.set_xticklabels(labels)
    ax.set_ylabel("information (bits)")
    save(fig, outdir, gene, name)
    rep.append((name, note, lines))


# ============================================================ 1. LENGTH =====
def plot_length(rows, dropped, outdir, gene, rep):
    n = len(rows)
    L = sorted(r["n"] for r in rows.values())
    null, g, p = comp_null(rows)
    med = int(np.median(L))
    cnt = Counter(L)

    fig, ax = plt.subplots(figsize=(4.4, 2.5))
    hi = min(max(L), 45)
    bins = np.arange(0, hi + 2) - 0.5
    ax.hist([x for x in L if x <= hi], bins=bins, color=C["hi"], edgecolor="none")
    if null:
        ax.axvline(null, color=C["null"], ls="--", lw=1.2)
        ax.text(null - 0.8, ax.get_ylim()[1] * .95,
                f"chance\n{null:.1f}", color=C["null"], fontsize=7,
                va="top", ha="right")
    ax.axvline(med, color=C["alt"], lw=1.2)
    ax.text(med + 0.8, ax.get_ylim()[1] * .95,
            f"observed\nmedian {med}", color=C["alt"], fontsize=7, va="top")
    ax.set_xlabel("extension length (codons)")
    ax.set_ylabel("species")
    ax.set_xlim(0, hi + 1)
    save(fig, outdir, gene, "length")

    # zoom on the bimodal gap
    fig, ax = plt.subplots(figsize=(3.2, 2.1))
    xs = list(range(8, 33))
    vals = [cnt[x] for x in xs]
    cols = [C["warn"] if v == 0 else C["hi"] for v in vals]
    ax.bar(xs, vals, color=cols, width=.8)
    ax.set_xlabel("extension length (codons)")
    ax.set_ylabel("species")
    ax.set_yscale("symlog", linthresh=1)
    save(fig, outdir, gene, "length_zoom")

    empty = [x for x in range(min(L), hi + 1) if cnt[x] == 0]
    rep.append(("length / length_zoom  [MAIN]", f"denominator A4 = {n}", [
        f"3'UTR GC (first 100 nt)     {100*g:.2f}%",
        f"P(stop) per random codon    {p:.6f}",
        f"chance expectation          {null:.2f} codons",
        f"observed median             {med} codons",
        f"ratio observed / chance     {med/null:.3f}",
        f"min {L[0]}  max {L[-1]}  mean {np.mean(L):.2f}",
        f"%% within 25-35 codons       {100*sum(1 for x in L if 25<=x<=35)/n:.2f}",
        f"species > {LONG_CUT} codons          {sum(1 for x in L if x>LONG_CUT)}",
        f"empty bins {min(L)}-{hi}: {empty if empty else 'none'}",
        f"frameshifted records removed: {dropped if dropped else 'none'}",
        "",
        "length : count : % of A4",
        *[f"  {k:>4} {cnt[k]:>5}  {100*cnt[k]/n:6.2f}%" for k in sorted(cnt)]]))


def plot_length_ranked(rows, outdir, gene, rep, top_n=10):
    """Ranked frequency of the commonest extension lengths.

    The histogram shows the shape of the distribution; this shows how much
    of the dataset a handful of discrete lengths accounts for.
    """
    n = len(rows)
    cnt = Counter(r["n"] for r in rows.values())
    top = cnt.most_common(top_n)
    covered = sum(c for _, c in top)

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    y = np.arange(len(top))[::-1]
    cols = []
    for L, _ in top:
        if L == 28:
            cols.append(C["hi"])          # canonical form
        elif L == 11:
            cols.append(C["alt"])         # truncated form
        else:
            cols.append(C["lo"])
    ax.barh(y, [100 * c / n for _, c in top], color=cols, height=.68)
    for yy, (L, c) in zip(y, top):
        ax.text(100 * c / n + 1.0, yy, f"{100*c/n:.1f}%  (n={c})",
                va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{L} codons" for L, _ in top], fontsize=7.5)
    ax.set_xlabel("species (%)")
    ax.set_xlim(0, max(100 * top[0][1] / n * 1.32, 25))
    save(fig, outdir, gene, "length_ranked")

    rep.append(("length_ranked  [MAIN]", f"denominator A4 = {n}", [
        f"the {top_n} commonest lengths account for {covered}/{n} "
        f"({100*covered/n:.2f}%) of species",
        f"distinct lengths observed: {len(cnt)}",
        "",
        f"{'length':>8}{'n':>6}{'%':>9}{'cumulative %':>15}",
        *[f"{L:>8}{c:>6}{100*c/n:>8.2f}%{100*sum(x for _,x in top[:i+1])/n:>14.2f}%"
          for i, (L, c) in enumerate(top)],
        "",
        "all observed lengths",
        *[f"  {k:>4} {cnt[k]:>5}  {100*cnt[k]/n:6.2f}%" for k in sorted(cnt)]]))


def plot_ecdf(rows, outdir, gene, rep):
    n = len(rows)
    L = np.sort([r["n"] for r in rows.values()])
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.step(L, np.arange(1, n + 1) / n, color=C["hi"], lw=1.4, where="post")
    null, _, _ = comp_null(rows)
    if null:
        x = np.arange(0, L.max() + 1)
        ax.plot(x, 1 - np.exp(-x / null), color=C["null"], ls="--", lw=1.1)
        ax.text(null * 2.3, .35, "geometric null", color=C["null"], fontsize=7)
    ax.set_xlabel("extension length (codons)")
    ax.set_ylabel("cumulative fraction")
    ax.set_xlim(0, min(L.max(), 45))
    save(fig, outdir, gene, "length_ecdf")
    rep.append(("length_ecdf  [SUPPLEMENTARY]", f"n = {n}",
                [f"median {int(np.median(L))}  q25 {int(np.percentile(L,25))}"
                 f"  q75 {int(np.percentile(L,75))}"]))


# ============================================================ 2. LOGOS ======
def plot_logos(rows, outdir, gene, rep, k=8):
    """Terminal and N-terminal logos, restricted to full-length extensions.

    The extension is length-variable, so a single anchoring would be
    misleading; sequences are aligned from the relevant terminus and only
    those long enough to fill the window are used.
    """
    seqs = [r["ext"] for r in rows.values() if len(r["ext"]) >= k]
    if len(seqs) < 20:
        rep.append(("logos", f"only {len(seqs)} sequences of >= {k} aa", []))
        return

    cterm = [s[-k:] for s in seqs]
    aa_logo(cterm, [str(-k + i) for i in range(k)], outdir, gene,
            "logo_cterm", rep,
            f"C-terminal {k} residues, aligned from the terminus; "
            f"n = {len(cterm)} sequences of >= {k} aa")

    nterm = [s[:k] for s in seqs]
    aa_logo(nterm, [str(i + 1) for i in range(k)], outdir, gene,
            "logo_nterm", rep,
            f"N-terminal {k} residues, aligned from the start; "
            f"n = {len(nterm)}")


# ======================================================== 3. MOTIF RANKING ==
def plot_motifs(rows, outdir, gene, rep):
    n = len(rows)
    hexa = Counter(r["ext"][-6:] for r in rows.values() if len(r["ext"]) >= 6)
    tri = Counter(r["ext"][-3:] for r in rows.values() if len(r["ext"]) >= 3)
    penta = Counter(r["ext"][:5] for r in rows.values() if len(r["ext"]) >= 5)

    def ranked(counter, name, colour, note, width=3.0):
        top = counter.most_common(6)
        fig, ax = plt.subplots(figsize=(width, 2.4))
        y = np.arange(len(top))[::-1]
        cols = [colour] + [C["lo"]] * (len(top) - 1)
        ax.barh(y, [100 * c / n for _, c in top], color=cols, height=.66)
        for yy, (_, c) in zip(y, top):
            ax.text(100 * c / n + 1.2, yy, f"{100*c/n:.1f}",
                    va="center", fontsize=7)
        ax.set_yticks(y)
        ax.set_yticklabels([m for m, _ in top], fontfamily="monospace", fontsize=7)
        ax.set_xlabel("species (%)")
        ax.set_xlim(0, max(100 * top[0][1] / n * 1.28, 20))
        save(fig, outdir, gene, name)
        rep.append((name, note, [f"{m:<10}{c:>5}  {100*c/n:6.2f}%"
                                 for m, c in counter.most_common(10)]))

    ranked(hexa, "motif_cterm", C["cterm"],
           f"C-terminal hexapeptide, denominator A4 = {n}  [MAIN]")
    ranked(penta, "motif_nterm", C["nterm"],
           f"N-terminal pentapeptide, denominator A4 = {n}  [SUPPLEMENTARY]")

    topC, tcC = hexa.most_common(1)[0]
    var = sum(v for m, v in hexa.items()
              if m != topC and len(m) == len(topC)
              and sum(a != b for a, b in zip(m, topC)) == 1)
    dli = sum(v for m, v in tri.items()
              if len(m) == 3 and m[0] == "D" and m[1] in "LI" and m[2] in "LI")
    rep.append(("motif summary", f"denominator A4 = {n}", [
        f"commonest hexapeptide {topC}: {tcC} ({100*tcC/n:.2f}%)",
        f"  + single-residue variants: {tcC+var} ({100*(tcC+var)/n:.2f}%)",
        f"D-[LI]-[LI] tripeptide: {dli} ({100*dli/n:.2f}%)",
        *[f"  tri {m:<6}{c:>5}  {100*c/n:6.2f}%" for m, c in tri.most_common(6)]]))
    return hexa, penta


# ========================================================= 4. ASYMMETRY =====
def plot_asymmetry(rows, hexa, penta, outdir, gene, rep):
    n = len(rows)
    topC, tcC = hexa.most_common(1)[0]
    topN, tcN = penta.most_common(1)[0]

    fig, ax = plt.subplots(figsize=(2.4, 2.4))
    v = [100 * tcN / n, 100 * tcC / n]
    ax.bar(["N-terminal\npentapeptide", "C-terminal\nhexapeptide"], v,
           color=[C["nterm"], C["cterm"]], width=.58)
    for i, (val, m) in enumerate(zip(v, (topN, topC))):
        ax.text(i, val + 2, f"{val:.1f}%", ha="center", fontsize=8, weight="bold")
        ax.text(i, val / 2, m, ha="center", fontsize=7,
                fontfamily="monospace", color="white")
    ax.set_ylabel("species carrying the commonest motif (%)")
    ax.set_ylim(0, 100)
    save(fig, outdir, gene, "asymmetry")
    rep.append(("asymmetry  [MAIN]", f"denominator A4 = {n}", [
        f"commonest C-terminal 6-mer {topC}: {100*tcC/n:.2f}%",
        f"commonest N-terminal 5-mer {topN}: {100*tcN/n:.2f}%",
        f"ratio: {tcC/tcN:.2f}x"]))


# ================================================= 5. PER-POSITION ==========
def plot_conservation(rows, outdir, gene, rep, modal=None):
    """Conservation at each position, full-length extensions only."""
    L = Counter(r["n"] for r in rows.values())
    modal = modal or L.most_common(1)[0][0]
    full = [r["ext"] for r in rows.values() if r["n"] == modal]
    n = len(full)
    if n < 20:
        rep.append(("conservation", f"only {n} at modal length {modal}", []))
        return

    frac, top = [], []
    for i in range(modal):
        col = Counter(s[i] for s in full)
        a, c = col.most_common(1)[0]
        frac.append(100 * c / n)
        top.append(a)

    fig, ax = plt.subplots(figsize=(5.2, 2.2))
    x = np.arange(modal)
    ax.bar(x, frac, color=C["hi"], width=.72)
    ax.set_xticks(x)
    ax.set_xticklabels(top, fontfamily="monospace", fontsize=6.5)
    ax.set_ylabel("commonest residue (%)")
    ax.set_xlabel(f"position within the {modal}-residue extension")
    ax.set_ylim(0, 105)
    save(fig, outdir, gene, "conservation")

    lines = [f"n = {n} sequences of exactly {modal} residues", "",
             f"{'pos':>4}{'aa':>5}{'%':>8}{'2nd':>6}{'%':>8}   S+T%"]
    for i in range(modal):
        col = Counter(s[i] for s in full)
        (a1, v1), *rest = col.most_common(2) + [("-", 0)]
        a2, v2 = rest[0]
        st = 100 * (col.get("S", 0) + col.get("T", 0)) / n
        flag = "  <- S/T interchangeable" if (a1 in "ST" and a2 in "ST"
                                              and v2 / v1 > 0.5) else ""
        lines.append(f"{i+1:>4}{a1:>5}{100*v1/n:>7.1f}{a2:>6}"
                     f"{100*v2/n:>7.1f}{st:>8.1f}{flag}")
    rep.append(("conservation  [SUPPLEMENTARY]",
                f"modal length {modal}, n = {n}", lines))


# ================================================ 6. CROSS-GENE NULL ========
def plot_null_vs_gc(genes, datadir, outdir, rep):
    """Chance expectation against 3'UTR GC, all genes — why a common
    threshold would be wrong."""
    pts = []
    for g in genes:
        try:
            rows, _ = a4(load(g, datadir))
        except SystemExit:
            continue
        except Exception:
            continue
        if not rows:
            continue
        null, gc, _ = comp_null(rows)
        med = int(np.median([r["n"] for r in rows.values()]))
        pts.append((g, 100 * gc, null, med, med / null, len(rows)))
    if len(pts) < 2:
        return

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for g, gc, null, med, ratio, n in pts:
        col = C["alt"] if g == "aqp4" else C["lo"]
        ax.scatter(gc, null, s=42, color=col, zorder=3)
        ax.scatter(gc, med, s=42, color=col, marker="^", zorder=3)
        ax.plot([gc, gc], [min(null, med), max(null, med)],
                color=col, lw=1.0, zorder=2)
        ax.annotate(g.upper(), (gc, max(null, med)), fontsize=7,
                    xytext=(0, 5), textcoords="offset points", ha="center")
    ax.set_xlabel("3'UTR GC content (%)")
    ax.set_ylabel("codons before first in-frame stop")
    ax.scatter([], [], s=42, color=C["lo"], label="chance expectation")
    ax.scatter([], [], s=42, color=C["lo"], marker="^", label="observed median")
    ax.legend(frameon=False, loc="upper left")
    save(fig, outdir, "all", "null_vs_gc")

    rep.append(("null_vs_gc  [SUPPLEMENTARY]",
                "chance expectation is gene-specific; a common threshold "
                "would misrank the genes", [
        f"{'gene':<7}{'n':>6}{'GC%':>8}{'null':>8}{'median':>8}{'ratio':>8}",
        *[f"{g:<7}{n:>6}{gc:>7.1f}{null:>8.1f}{med:>8}{ratio:>8.2f}"
          for g, gc, null, med, ratio, n in pts]]))


# ==================================================================== main ==
ALL = ["length", "length_zoom", "length_ranked", "length_ecdf",
       "logo_cterm", "logo_nterm", "motif_cterm", "motif_nterm",
       "asymmetry", "conservation", "null_vs_gc"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", default="aqp4")
    ap.add_argument("--datadir", default=".")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--only", default=None)
    ap.add_argument("--compare", nargs="*", default=None,
                    help="other genes, for the cross-gene null plot")
    ap.add_argument("--logo-window", type=int, default=8)
    args = ap.parse_args()

    gene = args.gene.lower()
    outdir = args.outdir or f"figures_{gene}_extension"
    os.makedirs(outdir, exist_ok=True)
    want = set(args.only.split(",")) if args.only else set(ALL)

    raw = load(gene, args.datadir)
    rows, dropped = a4(raw)
    print(f"=== {gene}: {len(raw)} with a 3'UTR record -> A4 = {len(rows)} "
          f"-> {outdir}/ ===")
    if dropped:
        print(f"  frameshifted, removed: {dropped}")
    rep = []

    if "length" in want or "length_zoom" in want:
        plot_length(rows, dropped, outdir, gene, rep)
    if "length_ranked" in want:
        plot_length_ranked(rows, outdir, gene, rep)
    if "length_ecdf" in want:
        plot_ecdf(rows, outdir, gene, rep)
    if "logo_cterm" in want or "logo_nterm" in want:
        plot_logos(rows, outdir, gene, rep, k=args.logo_window)

    hexa = penta = None
    if any(k in want for k in ("motif_cterm", "motif_nterm", "asymmetry")):
        hexa, penta = plot_motifs(rows, outdir, gene, rep)
    if "asymmetry" in want and hexa and penta:
        plot_asymmetry(rows, hexa, penta, outdir, gene, rep)
    if "conservation" in want:
        plot_conservation(rows, outdir, gene, rep)
    if "null_vs_gc" in want and args.compare is not None:
        plot_null_vs_gc([gene] + list(args.compare), args.datadir, outdir, rep)

    path = os.path.join(outdir, f"{gene}_extension_numbers.txt")
    with open(path, "w") as fh:
        fh.write(f"NUMBERS BEHIND EVERY PLOT — {gene.upper()} extension\n")
        fh.write("=" * 72 + "\n")
        fh.write(f"denominator A4 = {len(rows)} species: CDS + 3'UTR record,\n")
        fh.write("reading frame ending at a genuine in-frame stop codon,\n")
        fh.write("no frameshifted record.\n")
        for name, note, lines in rep:
            fh.write("\n" + "-" * 72 + f"\nPLOT: {name}\nNOTE: {note}\n\n")
            for l in lines:
                fh.write(f"  {l}\n")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
