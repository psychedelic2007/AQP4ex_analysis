#!/usr/bin/env python3
"""
Sequence logo for the AQP4 extension, without gap columns.

The problem
-----------
An MSA of a length-variable peptide contains many columns that are gaps in
most sequences. Plotting the alignment directly gives a logo dominated by
gap columns and positions that do not correspond to any real residue.

Three ways to get a continuous 1..28 logo, selectable with --mode:

  ref       Keep only the columns where a REFERENCE sequence has a residue.
            Gives exactly one position per reference residue, numbered in
            reference coordinates. Use this for a figure in human numbering.
            (default)

  modal     Ignore the alignment. Take every sequence whose ungapped length
            equals the modal length and stack them directly. No alignment
            assumptions at all, and the cleanest logo, but it discards
            sequences of other lengths.

  occupancy Keep alignment columns occupied by at least --min-occ of
            sequences. Useful for seeing whether any insertion column is
            well supported; positions are alignment columns, not residues.

Usage
-----
    python3 seqlogo.py --aln ext_aln.faa --mode ref --reference Homo_sapiens
    python3 seqlogo.py --aln ext_aln.faa --mode modal
    python3 seqlogo.py --fasta ext_prot.faa --mode modal        # unaligned input
    python3 seqlogo.py --aln ext_aln.faa --mode ref --offset 324 --out logo

--offset numbers the x-axis in protein coordinates: position i of the
reference is labelled offset + i. For the AQP4 extension use 324, so the
first residue is 325.

--pymol additionally writes a PyMOL script that colours a structure by the
per-position conservation computed here, storing the value in the B-factor
column so it can be recoloured or exported afterwards. Residue numbering in
the script follows --offset, so it must match the numbering of the PDB or
CIF file being coloured.

Gaps are never counted toward frequencies; each column is normalised over
the residues actually present, and the number of sequences contributing is
reported per position.
"""

import argparse
import os
import sys
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from matplotlib.font_manager import FontProperties

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("biopython required:  pip install biopython")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 150,
})

AA20 = "ACDEFGHIKLMNPQRSTVWY"

# Colour by physicochemical class. Override with --scheme mono for black.
SCHEMES = {
    "class": {**{a: "#1b7f79" for a in "AVLIMFWP"},      # hydrophobic
              **{a: "#5d708f" for a in "STNQCGY"},        # polar
              **{a: "#c53030" for a in "DE"},             # acidic
              **{a: "#2b6cb0" for a in "KRH"}},           # basic
    "mono": {a: "#222222" for a in AA20},
}


def letter(ch, x, y, w, h, ax, colour):
    """Draw one character scaled to fill (w, h) at (x, y)."""
    tp = TextPath((0, 0), ch, size=1,
                  prop=FontProperties(family="DejaVu Sans", weight="bold"))
    bb = tp.get_extents()
    if bb.width <= 0 or bb.height <= 0 or h <= 0:
        return
    t = (Affine2D().translate(-bb.x0, -bb.y0)
         .scale(w / bb.width, h / bb.height).translate(x, y))
    ax.add_patch(PathPatch(t.transform_path(tp), facecolor=colour,
                           edgecolor="none"))


# --------------------------------------------------------------- column sets
def cols_by_reference(recs, reference):
    """Alignment columns where the reference has a residue, in order."""
    ref = next((str(r.seq) for r in recs if r.id == reference), None)
    if ref is None:
        ids = [r.id for r in recs][:6]
        sys.exit(f"reference '{reference}' not in alignment. First ids: {ids}")
    cols = [i for i, c in enumerate(ref) if c not in "-."]
    return cols, "".join(ref[i] for i in cols)


def cols_by_occupancy(recs, min_occ):
    n = len(recs)
    L = len(recs[0].seq)
    cols = []
    for i in range(L):
        occ = sum(1 for r in recs if str(r.seq)[i] not in "-.") / n
        if occ >= min_occ:
            cols.append(i)
    return cols, None


def seqs_modal(recs):
    """Ungapped sequences of the modal length, stacked directly."""
    ungapped = [str(r.seq).replace("-", "").replace(".", "").upper()
                for r in recs]
    modal = Counter(len(s) for s in ungapped).most_common(1)[0][0]
    return [s for s in ungapped if len(s) == modal], modal


# ------------------------------------------------------------------ the logo
def draw(columns, labels, outdir, name, args, ref_seq=None):
    """columns: list of lists of residues (gaps already removed)."""
    L = len(columns)
    palette = SCHEMES[args.scheme]
    width = args.width or max(3.0, 0.24 * L + 0.9)

    fig, ax = plt.subplots(figsize=(width, args.height))
    peak, report = 0.0, []

    for i, col in enumerate(columns):
        col = [c for c in col if c in AA20]
        tot = len(col)
        if tot == 0:
            report.append((labels[i], 0, 0.0, []))
            continue
        cnt = Counter(col)
        freqs = {a: cnt.get(a, 0) / tot for a in AA20}
        H = -sum(p * np.log2(p) for p in freqs.values() if p > 0)
        if args.scale == "bits":
            R = max(np.log2(20) - H - (19 / (2 * np.log(2) * tot)), 0.0)
            ymax = np.log2(20)
        else:                                   # probability stack, height 1
            R, ymax = 1.0, 1.0
        peak = max(peak, R)

        y = 0.0
        for a, p in sorted(freqs.items(), key=lambda t: t[1]):
            if p < args.min_freq:
                continue
            h = p * R
            letter(a, i + 0.10, y, 0.80, h, ax, palette.get(a, "#666"))
            y += h
        report.append((labels[i], tot, R, cnt.most_common(4)))

    ax.set_xlim(0, L)
    ax.set_ylim(0, (min(peak * 1.12, ymax) if args.scale == "bits" else 1.02))
    ax.set_xticks(np.arange(L) + 0.5)
    ax.set_xticklabels(labels, rotation=90 if args.rotate else 0, fontsize=6.5)
    ax.set_ylabel("information (bits)" if args.scale == "bits"
                  else "frequency")
    if ref_seq and args.show_reference:
        for i, c in enumerate(ref_seq):
            ax.text(i + 0.5, -0.06 * ax.get_ylim()[1], c, ha="center",
                    va="top", fontsize=6.5, fontfamily="monospace",
                    color="#666", transform=ax.transData)
    for ext in ("svg", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"),
                    bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {name}.svg / .pdf")
    return report


# ------------------------------------------------------------------ pymol --
def write_pymol(report, path, args, mode):
    """Emit a PyMOL script colouring a structure by per-position conservation.

    Values are written into the B-factor column, which is the standard route:
    it survives saving, can be recoloured with any spectrum, and lets the
    numbers be recovered from the modified structure.

    Two scores are offered because they answer different questions:
      bits     information content, 0 to log2(20); rewards a column dominated
               by one residue regardless of which
      identity frequency of the commonest residue, 0-100; easier to explain
               and directly comparable to a percent-identity figure
    """
    rows = [(lab, tot, R, top) for lab, tot, R, top in report if tot > 0]
    if not rows:
        print("  [pymol] no positions with data — script not written")
        return

    if args.pymol_score == "bits":
        vals = {int(lab): R for lab, tot, R, top in rows}
        lo, hi = 0.0, float(np.log2(20))
        label = "information content (bits)"
    else:
        vals = {int(lab): 100 * top[0][1] / tot for lab, tot, R, top in rows}
        lo, hi = 0.0, 100.0
        label = "identity of the commonest residue (%)"

    vmin, vmax = min(vals.values()), max(vals.values())
    chain = args.pymol_chain
    sel = f"chain {chain} and " if chain else ""
    obj = args.pymol_object or "structure"

    with open(path, "w") as fh:
        fh.write(f"""# PyMOL script — conservation of the {args.out} positions
# generated by seqlogo.py (mode={mode}, score={args.pymol_score})
#
# Score: {label}
# Range in this dataset: {vmin:.3f} to {vmax:.3f}
# Positions: {min(vals)}-{max(vals)} ({len(vals)} residues)
#
# USAGE
#   1. Load your structure first, or edit the `load` line below.
#   2. @{os.path.basename(path)}
#
# IMPORTANT
#   Residue numbering below follows --offset {args.offset}. It must match the
#   numbering in your structure file. If your model numbers the extension
#   from 1, rerun with --offset 0, or renumber the structure:
#       alter {sel}polymer, resi=str(int(resi)+{args.offset})
#   Positions absent from the alignment are left at B = -1 and coloured grey,
#   so unmodelled or unscored residues are visibly distinct from conserved ones.

# load {obj}.pdb, {obj}          # uncomment and edit if not already loaded

# --- reset B-factors so unscored residues are identifiable ----------------
alter {obj}, b = -1.0

# --- per-residue conservation --------------------------------------------
""")
        for r in sorted(vals):
            fh.write(f"alter {obj} and {sel}resi {r}, b = {vals[r]:.4f}\n")

        fh.write(f"""
# --- display --------------------------------------------------------------
hide everything, {obj}
show cartoon, {obj}
set cartoon_transparency, 0.0

# grey for anything without a score
color grey70, {obj} and b < 0

# spectrum over the scored range only
spectrum b, blue_white_red, {obj} and b > -0.5, minimum={lo}, maximum={hi}

# highlight the scored region
select scored, {obj} and {sel}resi {min(vals)}-{max(vals)}
show sticks, scored and not (name C+N+O)
set cartoon_side_chain_helper, 1
deselect

bg_color white
set ray_opaque_background, 0
orient scored

# --- optional: label the most conserved positions -------------------------
# select topcons, {obj} and b > {vmin + 0.75*(vmax-vmin):.3f}
# label topcons and name CA, "%s%s" % (one_letter[resn], resi)
# set label_size, 14

# --- optional: write the B-factors back out -------------------------------
# save {obj}_conservation.pdb, {obj}

print "conservation loaded: {label}"
print "scored residues {min(vals)}-{max(vals)}, B range {vmin:.3f}-{vmax:.3f}"
print "unscored residues have B = -1 and are coloured grey70"
""")
    print(f"  wrote {path}")
    print(f"  [pymol] score={args.pymol_score}  residues {min(vals)}-{max(vals)}"
          f"  B range {vmin:.3f}-{vmax:.3f}")

    tsv = path.rsplit(".", 1)[0] + "_bfactors.tsv"
    with open(tsv, "w") as fh:
        fh.write("residue\tn\tbits\tidentity_pct\ttop_residue\n")
        for lab, tot, R, top in rows:
            fh.write(f"{lab}\t{tot}\t{R:.4f}\t{100*top[0][1]/tot:.2f}\t"
                     f"{top[0][0]}\n")
    print(f"  wrote {tsv}")


# ------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--aln", help="aligned FASTA")
    src.add_argument("--fasta", help="unaligned FASTA (forces --mode modal)")
    ap.add_argument("--mode", default="ref", choices=["ref", "modal", "occupancy"])
    ap.add_argument("--reference", default="Homo_sapiens")
    ap.add_argument("--offset", type=int, default=0,
                    help="label position i as offset+i (324 for the AQP4 extension)")
    ap.add_argument("--min-occ", type=float, default=0.5,
                    help="occupancy mode: minimum fraction of non-gap residues")
    ap.add_argument("--min-freq", type=float, default=0.01,
                    help="do not draw letters below this frequency")
    ap.add_argument("--scale", default="bits", choices=["bits", "prob"])
    ap.add_argument("--scheme", default="class", choices=["class", "mono"])
    ap.add_argument("--width", type=float, default=None)
    ap.add_argument("--height", type=float, default=2.2)
    ap.add_argument("--rotate", action="store_true", help="rotate x labels")
    ap.add_argument("--show-reference", action="store_true",
                    help="print the reference sequence beneath the axis")
    ap.add_argument("--out", default="seqlogo")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--pymol", action="store_true",
                    help="also write a PyMOL script colouring a structure by "
                         "per-position conservation")
    ap.add_argument("--pymol-score", default="identity",
                    choices=["identity", "bits"],
                    help="quantity written to the B-factor column")
    ap.add_argument("--pymol-object", default=None,
                    help="PyMOL object name (default: structure)")
    ap.add_argument("--pymol-chain", default=None,
                    help="restrict colouring to one chain")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    path = args.aln or args.fasta
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs:
        sys.exit(f"no sequences in {path}")
    mode = "modal" if args.fasta else args.mode

    if mode == "modal":
        seqs, modal = seqs_modal(recs)
        if len(seqs) < 10:
            sys.exit(f"only {len(seqs)} sequences at modal length {modal}")
        columns = [[s[i] for s in seqs] for i in range(modal)]
        labels = [str(args.offset + i + 1) for i in range(modal)]
        ref_seq = None
        print(f"{len(recs)} sequences -> {len(seqs)} of modal length {modal}")
    else:
        L = {len(r.seq) for r in recs}
        if len(L) != 1:
            sys.exit(f"--mode {mode} needs an alignment; sequence lengths vary {sorted(L)[:5]}")
        if mode == "ref":
            cols, ref_seq = cols_by_reference(recs, args.reference)
            labels = [str(args.offset + i + 1) for i in range(len(cols))]
            print(f"{len(recs)} sequences, {list(L)[0]} columns -> "
                  f"{len(cols)} columns where {args.reference} has a residue")
        else:
            cols, ref_seq = cols_by_occupancy(recs, args.min_occ)
            labels = [str(c + 1) for c in cols]
            print(f"{len(recs)} sequences, {list(L)[0]} columns -> "
                  f"{len(cols)} columns with occupancy >= {args.min_occ}")
        columns = [[str(r.seq)[c].upper() for r in recs] for c in cols]

    report = draw(columns, labels, args.outdir, args.out, args, ref_seq)

    txt = os.path.join(args.outdir, f"{args.out}_numbers.txt")
    with open(txt, "w") as fh:
        fh.write(f"SEQUENCE LOGO — {os.path.basename(path)}\n")
        fh.write(f"mode={mode}  scale={args.scale}  "
                 f"sequences={len(recs)}  positions={len(columns)}\n")
        if mode == "ref":
            fh.write(f"reference={args.reference}\n")
            fh.write("Columns gapped in the reference are excluded, so each\n")
            fh.write("position corresponds to one reference residue.\n")
        elif mode == "modal":
            fh.write("Alignment not used; sequences of modal ungapped length\n")
            fh.write("stacked directly.\n")
        fh.write("Gaps are excluded from frequencies; n is the number of\n")
        fh.write("residues actually present at each position.\n")
        fh.write("\n" + "-" * 66 + "\n")
        fh.write(f"{'pos':>6}{'n':>7}{'bits':>8}   composition\n")
        for lab, tot, R, top in report:
            comp = "  ".join(f"{a}={100*c/tot:.1f}%" for a, c in top) if tot else "-"
            fh.write(f"{lab:>6}{tot:>7}{R:>8.3f}   {comp}\n")
    print(f"  wrote {txt}")

    if args.pymol:
        write_pymol(report, os.path.join(args.outdir, f"{args.out}.pml"),
                    args, mode)

    if mode == "ref":
        print("\nNote: 'ref' mode shows only positions present in the reference.")
        print("Insertions carried by other species are not displayed. Use")
        print("--mode occupancy if you need to see whether any are well supported.")


if __name__ == "__main__":
    main()
