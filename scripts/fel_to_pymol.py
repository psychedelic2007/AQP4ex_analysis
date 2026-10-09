#!/usr/bin/env python3
"""
PyMOL colouring scripts for site-wise FEL results, per gene and region.

    python3 fel_to_pymol.py --gene aqp4 --region tail
    python3 fel_to_pymol.py --gene aqp1 --region channel --pdb alphafold/aqp1_channel.pdb
    python3 fel_to_pymol.py --gene all  --region both          # paralogue panel set
    python3 fel_to_pymol.py --gene aqp4 --region both --quantity alpha   # or beta
    python3 fel_to_pymol.py --gene all  --region both --quantity all     # omega, alpha, beta

Source of the numbers
---------------------
Site tables are rebuilt from the joint channel+extension FEL fit of record
(results_final/<gene>/fel/joint_mafft.json) with final_analysis.analyse_job, so
the occupancy rule, alpha limits, BH q and calls are identical to
results_final/final_results.txt. For AQP4 the rebuilt tables reproduce
results_final/tables/aqp4_joint_mafft_{channel,extension}.tsv exactly. Paralogue
channel tables are not written by final_analysis.py; this script writes them to
the output directory. --table overrides everything with an explicit TSV.

What gets coloured
------------------
Only TESTED sites (call != "not tested") are coloured. Untested sites (low
occupancy, alpha > 10, alpha at the 0 boundary) are grey: their omega is
undefined or unreliable, e.g. 92755 at AQP4 L325 in the tetrapod subset.
Alignment columns that are gaps in the human sequence (ins>N) have no residue
in a human structure and are dropped. For AQP2, AQP6 and AQP8 some of these
are tested tail sites, so the tail figure shows fewer tested sites than
final_results.txt reports. The script prints them; name them in the legend.

Scales: the same for every gene and region
------------------------------------------
Paralogue panels only compare if every panel uses the same mapping, so no scale
is fitted to the data being drawn.
  omega       log2(omega), clipped to +/- --omega-cap (default 4, i.e. 1/16-16);
              red = purifying, white = 1 (neutral), blue = diversifying.
              omega is a ratio, and a linear 0-2 ramp would squeeze all
              purifying sites into half the colour range. The cap of 4 keeps
              the channel (median omega ~0.1) from saturating to flat red.
  alpha, beta one shared linear scale 0..--rate-max. By default this is the
              95th percentile of tested alpha and beta pooled over all six
              genes and both regions, so it does not change with --gene.
              Values above the maximum are clipped and counted in the header.
The B-factor column holds the DISPLAYED value (log2 omega, or the clipped
rate). The raw values are in the *_sites.tsv file written next to each script.
colorbar_<quantity>.png (horizontal, needs matplotlib) is drawn from the same
palette and limits as the spectrum command; the omega bar is in log2 units, and
its arrow ends mean values beyond the cap are clipped to the end colour.

Colour is not the test
----------------------
A tested but non-significant site can still have an extreme point-estimate
omega when alpha is small. Significant sites (q <= 0.20, the manuscript
threshold; --bonferroni for the stricter set) are drawn as sticks. Say in the
legend that colour is the FEL point estimate and sticks mark significance.

Numbering
---------
Channel residues use human reference numbering (1..N). Extension residues
continue past the stop codon (AQP4: first extension residue 325), so a tail
map needs a read-through model with that numbering. Pass --pdb and the script
checks that the residue at each number matches the human amino acid in the
table, and refuses to write if fewer than 90% match.
"""

import argparse
import csv
import math
import os
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/
ROOT = os.path.dirname(HERE)                                 # repository root
GENES = ("aqp1", "aqp2", "aqp4", "aqp5", "aqp6", "aqp8")
REGION_ALIAS = {"tail": "extension", "extension": "extension", "channel": "channel"}
SIG_CALLS = ("purifying", "diversifying")

SPEC = {
    "omega": dict(
        palette="red_white_blue",
        meaning="beta/alpha, shown as log2",
        reading="red = purifying (omega < 1), white = neutral, blue = diversifying"),
    "alpha": dict(
        palette="white_blue",
        meaning="synonymous substitution rate",
        reading="blue = substitutions reached this site freely. Alpha is the "
                "internal control: low beta is evidence only where alpha is high"),
    "beta": dict(
        palette="red_white",
        meaning="non-synonymous substitution rate",
        reading="red = amino-acid changes were removed"),
}

AA3 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G",
           HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S",
           THR="T", TRP="W", TYR="Y", VAL="V", SEC="U", MSE="M", HID="H", HIE="H",
           HIP="H", CYX="C")


# ------------------------------------------------------------------ loading --
def _fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def rows_from_tsv(path):
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        need = {"human_residue", "human_aa", "alpha", "beta", "omega", "q", "call"}
        if not need <= set(rdr.fieldnames or []):
            sys.exit(f"{path}: need columns {sorted(need)}; found {rdr.fieldnames}")
        return [dict(label=r["human_residue"], human=r["human_aa"],
                     alpha=_fnum(r["alpha"]), beta=_fnum(r["beta"]),
                     omega=_fnum(r["omega"]), q=_fnum(r["q"]),
                     bonf=r.get("bonferroni", "0") in ("1", "True"),
                     call=r["call"], reason=r.get("exclusion_reason", ""))
                for r in rdr]


def rows_from_fit(gene, region):
    """Rebuild one module's site table from the joint fit of record."""
    cwd = os.getcwd()
    os.chdir(ROOT)                       # final_analysis uses paths relative to the repo root
    sys.path.insert(0, HERE)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import final_analysis as fa
        d = f"{fa.OUT}/{gene}"
        js, aln = f"{d}/fel/joint_mafft.json", f"{d}/fel/joint_mafft.fna"
        for p in (js, aln, f"{d}/channel_codon.fna"):
            if not os.path.exists(p):
                sys.exit(f"missing {os.path.join(ROOT, p)}; run scripts/final_analysis.py first")
        ch = fa.read_fasta(f"{d}/channel_codon.fna")
        nchan = len(next(iter(ch.values()))) // 3
        mod = fa.analyse_job(dict(gene=gene, json=js, aln=aln, nchan=nchan,
                                  module="joint"))["modules"][region]
    finally:
        os.chdir(cwd)
    return [dict(label=r["label"], human=r["human"], alpha=r["alpha"], beta=r["beta"],
                 omega=r["omega"], q=r["q"], bonf=bool(r["bonf"]), call=r["call"],
                 reason=r["reason"]) for r in mod]


def residues(rows, shift=0):
    """Keep rows that map to a human residue; return (kept, dropped insertion rows).

    ref = human reference number (as in the tables); res = structure number
    (ref + shift). Everything written to PyMOL uses res.
    """
    kept, ins = [], []
    for r in rows:
        if r["label"].isdigit():
            kept.append(dict(r, ref=int(r["label"]), res=int(r["label"]) + shift))
        else:
            ins.append(r)
    return kept, ins


def tested(rows):
    return [r for r in rows if r["call"] != "not tested"]


# ----------------------------------------------------------------- mapping --
def display_value(r, quantity, omega_cap, rate_max):
    if quantity == "omega":
        w = r["omega"]
        if not w == w:                                   # nan
            return None, False
        if w <= 0:
            return -omega_cap, True
        v = math.log2(w)
        return max(-omega_cap, min(omega_cap, v)), abs(v) > omega_cap
    v = r[quantity]
    if not v == v:
        return None, False
    return min(v, rate_max), v > rate_max


def pooled_rate_max(pct=95):
    vals = []
    for g in GENES:
        for reg in ("channel", "extension"):
            try:
                for r in tested(rows_from_fit(g, reg)):
                    vals += [r["alpha"], r["beta"]]
            except SystemExit:
                print(f"  [warn] {g} {reg}: fit unavailable, left out of the rate scale")
    if not vals:
        sys.exit("no tested sites anywhere; cannot set a rate scale")
    vals.sort()
    k = (len(vals) - 1) * pct / 100
    lo, hi = int(math.floor(k)), int(math.ceil(k))
    return vals[lo] + (vals[hi] - vals[lo]) * (k - lo), len(vals)


def compress(nums):
    """[1,2,3,7,9,10] -> '1-3+7+9-10' (PyMOL resi syntax)."""
    nums, out, i = sorted(set(nums)), [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if i == j else f"{nums[i]}-{nums[j]}")
        i = j + 1
    return "+".join(out)


# --------------------------------------------------------------- PDB check --
def pdb_sequence(path, chain):
    seq = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
                if chain and line[21] != chain:
                    continue
                try:
                    resi = int(line[22:26])
                except ValueError:
                    continue
                seq.setdefault(resi, AA3.get(line[17:20].strip(), "X"))
    return seq


def check_numbering(rows, pdb, chain, force):
    if pdb.lower().endswith((".cif", ".mmcif")):
        print("  [warn] numbering check reads PDB format only; convert or check by hand")
        return
    seq = pdb_sequence(pdb, chain)
    if not seq:
        sys.exit(f"{pdb}: no CA atoms" + (f" in chain {chain}" if chain else ""))
    present = [r for r in rows if r["res"] in seq]
    match = [r for r in present if seq[r["res"]] == r["human"]]
    missing = len(rows) - len(present)
    frac = len(match) / len(present) if present else 0.0
    print(f"  numbering vs {os.path.basename(pdb)}: {len(match)}/{len(present)} residues "
          f"match the human amino acid ({100 * frac:.0f}%), {missing} not in structure")
    bad = [r for r in present if seq[r["res"]] != r["human"]][:8]
    for r in bad:
        print(f"     resi {r['res']}: table {r['human']}, structure {seq[r['res']]}")
    if (not present or frac < 0.90) and not force:
        best = best_shift(rows, seq)
        if best:
            k, n, m = best
            cur = rows[0]["res"] - rows[0]["ref"]
            sys.exit(f"  numbering does not match. Shift {cur + k:+d} gives {n}/{m} matches: "
                     f"rerun with --resi-shift {cur + k} (table {rows[0]['ref']} -> "
                     f"structure {rows[0]['ref'] + cur + k})")
        sys.exit("  numbering does not match; renumber the structure or pass --force")


def best_shift(rows, seq):
    """Find the offset (structure - current numbering) that best matches the sequence."""
    best = None
    for k in range(min(seq) - max(r["res"] for r in rows),
                   max(seq) - min(r["res"] for r in rows) + 1):
        pres = [r for r in rows if r["res"] + k in seq]
        if len(pres) < max(5, len(rows) // 2):
            continue
        n = sum(seq[r["res"] + k] == r["human"] for r in pres)
        if n / len(pres) >= 0.90 and (best is None or n > best[1]):
            best = (k, n, len(pres))
    return best


# ---------------------------------------------------------------- writing --
def write_script(rows, gene, region, quantity, args, outdir, scale, ins, src):
    spec = SPEC[quantity]
    obj = args.object or gene
    sel = f" and chain {args.chain}" if args.chain else ""
    lo, hi = scale
    shift = rows[0]["res"] - rows[0]["ref"]
    pool = tested(rows)
    if args.sig_only:
        pool = [r for r in pool if is_sig(r, args)]
    vals, clipped = {}, 0
    for r in pool:
        v, c = display_value(r, quantity, args.omega_cap, args.rate_max)
        if v is None:
            continue
        vals[r["res"]] = v
        clipped += c
    if not vals:
        print(f"  [skip] {gene} {region} {quantity}: nothing to colour")
        return None
    sig = sorted(r["res"] for r in tested(rows) if is_sig(r, args))
    sticks = args.sticks == "on" or (args.sticks == "auto" and region == "extension")
    raw = [r[quantity] for r in pool if r[quantity] == r[quantity]]
    tag = f"{gene}_{'tail' if region == 'extension' else 'channel'}_{quantity}"
    path = os.path.join(outdir, f"{tag}.pml")
    sig_rule = "Bonferroni" if args.bonferroni else "BH q <= 0.20"
    ramp = ("log2(omega) from {:+.1f} to {:+.1f} (omega {:.3g} to {:.3g}), white at "
            "omega = 1".format(lo, hi, 2 ** lo, 2 ** hi) if quantity == "omega"
            else f"{lo:.3f} to {hi:.3f}, shared by alpha and beta")

    fh = _Pml()
    if True:
        fh.write(f"""# PyMOL: {gene.upper()} {region}, coloured by {quantity} ({spec['meaning']})
# generated by fel_to_pymol.py from {src}
#
# Reading: {spec['reading']}
# Colour scale : {ramp}
#                identical for every gene and region, so panels compare directly
# Raw range    : {min(raw):.3g} to {max(raw):.3g} over {len(vals)} coloured sites;
#                {clipped} clipped at the scale limit
# Rows         : {len(rows)} human residues ({len(tested(rows))} tested). {ins[0]} columns
#                are gaps in human and are not drawn; {ins[1]} of them were tested
#                and {ins[2]} significant, so they are missing from this figure.
# Grey         : not tested{' or not significant' if args.sig_only else ''}
# Sticks       : {f"{len(sig)} significant sites ({sig_rule})" if sticks else f"off ({len(sig)} significant, {sig_rule}; --sticks on to draw)"}
# Numbering    : structure resi = human reference number {shift:+d}
#                ({rows[0]['ref']} -> {rows[0]['res']}), set by --resi-shift
# B-factor     : the displayed value, not the raw rate; raw values in {tag}_sites.tsv
#
# USAGE (in PyMOL)
#   load <structure>.pdb, {obj}     # skip if --pdb was given
#   @{os.path.basename(path)}
""")
        if args.pdb:
            fh.write(f"\ndelete {obj}\nload {os.path.abspath(args.pdb)}, {obj}\n")
        fh.write(f"\nalter {obj}, b = 0.0\n")
        for r in sorted(vals):
            fh.write(f"alter {obj}{sel} and resi {r}, b = {vals[r]:.4f}\n")
        fh.write(f"""
rebuild
hide everything, {obj}
show cartoon, {obj}
bg_color white
set ray_opaque_background, 0
set cartoon_side_chain_helper, 1

select fel_scored, {obj}{sel} and resi {compress(vals)}
color grey70, {obj}
spectrum b, {spec['palette']}, fel_scored, minimum={lo}, maximum={hi}
""")
        if sig:
            fh.write(f"""
select fel_significant, {obj}{sel} and resi {compress(sig)}
{"" if sticks else "# "}show sticks, fel_significant and not name N+C+O
""")
        fh.write(f"""deselect
orient {obj}{sel} and resi {min(vals)}-{max(vals)}

{pml_print(f"{gene.upper()} {region} {quantity}: {ramp}")}
{pml_print(f"{len(vals)} coloured / {clipped} clipped / {len(sig)} significant ({sig_rule})")}

# --- optional -------------------------------------------------------------
# ray 1600, 1200
# png {tag}.png, dpi=300
""")
    fh.save(path)

    with open(os.path.join(outdir, f"{tag}_sites.tsv"), "w") as fh:
        fh.write("structure_resi\thuman_residue\thuman_aa\talpha\tbeta\tomega\tq\tbonferroni\tcall\t"
                 "displayed_value\texclusion_reason\n")
        for r in rows:
            v = vals.get(r["res"])
            fh.write(f"{r['res']}\t{r['ref']}\t{r['human']}\t{r['alpha']:.5g}\t{r['beta']:.5g}\t"
                     f"{r['omega']:.5g}\t{r['q']:.4g}\t{int(r['bonf'])}\t{r['call']}\t"
                     f"{'' if v is None else f'{v:.4f}'}\t{r['reason']}\n")
    print(f"  wrote {os.path.basename(path)}: {len(vals)} coloured, {clipped} clipped, "
          f"{len(sig)} significant")
    return path


# PyMOL spectrum palettes are linear RGB ramps through its named colours
# (red = 1,0,0 / white = 1,1,1 / blue = 0,0,1); the colour bar uses the same anchors.
PALETTE_RGB = {
    "red_white_blue": [(1, 0, 0), (1, 1, 1), (0, 0, 1)],
    "white_blue": [(1, 1, 1), (0, 0, 1)],
    "red_white": [(1, 0, 0), (1, 1, 1)],
}


def write_colorbar(quantity, scale, outdir, labelsize=30):
    """Horizontal colour bar with the palette and limits given to PyMOL's spectrum.

    Arrow ends mark that values beyond the limits are clipped to the end colour.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.colorbar as mcolorbar
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [warn] matplotlib not installed; no colour bar written")
        return None
    lo, hi = scale
    if not math.isfinite(hi):
        return None
    cmap = mcolors.LinearSegmentedColormap.from_list(
        SPEC[quantity]["palette"], PALETTE_RGB[SPEC[quantity]["palette"]])
    fig, ax = plt.subplots(figsize=(6, 1.5))
    cb = mcolorbar.ColorbarBase(ax, cmap=cmap, norm=mcolors.Normalize(vmin=lo, vmax=hi),
                                orientation="horizontal",
                                extend="both" if quantity == "omega" else "max")
    if quantity == "omega":
        cb.set_ticks([t for t in range(math.ceil(lo), math.floor(hi) + 1)
                      if t % 2 == 0])                       # log2(omega); 0 = neutral
    else:                                                   # label the (unrounded) maximum
        cb.set_ticks([lo, (lo + hi) / 2, hi])
        cb.set_ticklabels([f"{t:.2g}" for t in (lo, (lo + hi) / 2, hi)])
    cb.ax.tick_params(labelsize=labelsize)
    path = os.path.join(outdir, f"colorbar_{quantity}.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {os.path.basename(path)}: {SPEC[quantity]['palette']}, {lo:g} to {hi:g}")
    return path


class _Pml:
    """Collects the script text; semicolons are PyMOL command separators, even
    inside # comments, so they are replaced before the file is written."""

    def __init__(self):
        self.parts = []

    def write(self, text):
        self.parts.append(text)

    def save(self, path):
        with open(path, "w") as out:
            out.write("".join(self.parts).replace(";", ","))


def pml_print(text):
    """PyMOL's print evaluates comma-separated arguments; keep it to one string."""
    return 'print "' + text.replace(",", " -").replace('"', "'") + '"'


def is_sig(r, args):
    if r["call"] not in SIG_CALLS:
        return False
    return (r["bonf"] and r["beta"] < r["alpha"]) if args.bonferroni else True


# -------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True, choices=GENES + ("all",))
    ap.add_argument("--region", required=True,
                    choices=("tail", "extension", "channel", "both"))
    ap.add_argument("--quantity", default="omega", choices=("omega", "alpha", "beta", "all"),
                    help="default omega; 'all' writes omega, alpha and beta")
    ap.add_argument("--table", help="explicit site TSV (one gene and region only)")
    ap.add_argument("--pdb", help="structure: checks numbering and is loaded by the script")
    ap.add_argument("--resi-shift", type=int, default=0,
                    help="structure resi = table residue + SHIFT. A tail-only model "
                         "numbered from 1 needs -324 for AQP4 (first tail residue 325)")
    ap.add_argument("--force", action="store_true", help="write despite a numbering mismatch")
    ap.add_argument("--object", help="PyMOL object name (default: the gene name)")
    ap.add_argument("--chain", help="colour only this chain (default: all chains)")
    ap.add_argument("--outdir", default=os.path.join(ROOT, "figures", "pymol"))
    ap.add_argument("--sig-only", action="store_true",
                    help="colour only significant sites; the rest stay grey")
    ap.add_argument("--bonferroni", action="store_true",
                    help="significant = Bonferroni purifying instead of q <= 0.20")
    ap.add_argument("--sticks", default="auto", choices=("auto", "on", "off"),
                    help="significant sites as sticks; auto = tail only (in the channel "
                         "~95%% of sites are significant, so sticks would cover it)")
    ap.add_argument("--omega-cap", type=float, default=4.0,
                    help="clip |log2 omega| here (default 4 = omega 0.0625-16)")
    ap.add_argument("--rate-max", type=float,
                    help="alpha/beta scale maximum (default: pooled 95th percentile)")
    args = ap.parse_args()

    genes = GENES if args.gene == "all" else (args.gene,)
    regions = (("channel", "extension") if args.region == "both"
               else (REGION_ALIAS[args.region],))
    quants = ("omega", "alpha", "beta") if args.quantity == "all" else (args.quantity,)
    if args.table and (len(genes) > 1 or len(regions) > 1):
        sys.exit("--table takes one gene and one region")
    if args.pdb and len(genes) > 1:
        sys.exit("--pdb is one structure; give one --gene with it")
    if args.pdb and not os.path.exists(args.pdb):
        sys.exit(f"missing: {args.pdb}")

    if {"alpha", "beta"} & set(quants) and args.rate_max is None:
        args.rate_max, n = pooled_rate_max()
        print(f"alpha/beta scale 0-{args.rate_max:.3f} (95th percentile of {n} tested "
              f"alpha and beta values, all genes and regions)")
    if args.rate_max is None:
        args.rate_max = float("inf")

    os.makedirs(args.outdir, exist_ok=True)
    written = []
    for g in genes:
        for reg in regions:
            if args.table:
                rows, src = rows_from_tsv(args.table), os.path.basename(args.table)
            else:
                rows = rows_from_fit(g, reg)
                src = f"results_final/{g}/fel/joint_mafft.json ({reg} module)"
            rows, dropped = residues(rows, args.resi_shift)
            t = tested(rows)
            npur = sum(r["call"] == "purifying" for r in t)
            ndiv = sum(r["call"] == "diversifying" for r in t)
            print(f"{g} {reg}: {len(rows)} residues, {len(t)} tested, {npur} purifying, "
                  f"{ndiv} diversifying (q <= 0.20); {len(dropped)} insertion columns dropped")
            lost = tested(dropped)
            if lost:
                ls = [r for r in lost if r["call"] in SIG_CALLS]
                print(f"  [NOTE] {len(lost)} tested columns ({len(ls)} significant: "
                      f"{', '.join(r['label'] + '/' + r['call'][:3] for r in ls) or 'none'}) "
                      f"are gaps in human and cannot be shown on a human structure. "
                      f"Report them in the legend.")
            ins = (len(dropped), len(lost), sum(r["call"] in SIG_CALLS for r in lost))
            if args.pdb:
                check_numbering(rows, args.pdb, args.chain, args.force)
            for q in quants:
                scale = ((-args.omega_cap, args.omega_cap) if q == "omega"
                         else (0.0, args.rate_max))
                p = write_script(rows, g, reg, q, args, args.outdir, scale, ins, src)
                if p:
                    written.append(p)

    # one bar per quantity: the scale is the same for every gene and region
    for q in quants:
        write_colorbar(q, (-args.omega_cap, args.omega_cap) if q == "omega"
                       else (0.0, args.rate_max), args.outdir)

    print(f"\n{len(written)} scripts in {args.outdir}")
    print("""
IN PYMOL
  load <structure>.pdb, <gene>        # object name must match (see --object)
  @<gene>_<region>_<quantity>.pml

FOR THE LEGEND
  Colour is the FEL point estimate on one fixed scale for every gene and region;
  sticks mark significant sites. Grey = not tested (occupancy < 50%, alpha > 10,
  or alpha at the boundary), where omega is undefined or unreliable.
  Sites where beta ~ 0 have omega ~ 0 and sit at the red end of the scale
  (counted as clipped); that is the expected state of a conserved channel.""")


if __name__ == "__main__":
    main()
