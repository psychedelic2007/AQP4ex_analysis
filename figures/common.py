"""Shared settings, palette and data loaders for the AQP4ex manuscript figures.

Every panel script in figures/panels/ imports this module. All data are read
from the analysis of record (results_final/, produced by scripts/final_analysis.py) or
from the raw NCBI downloads in data/; nothing is hard-coded except
values that final_analysis.py only writes to results_final/final_results.txt
(the Mk-model fits), which are parsed from that report.
"""
import json
import math
import os
import re
import warnings
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)

# ------------------------------------------------------------------ paths --
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results_final")
DATA = os.path.join(ROOT, "data")
SCRIPTS = os.path.join(ROOT, "scripts")      # final_analysis.py is imported from here
TAB = os.path.join(RES, "tables")
OUTDIR = os.path.join(ROOT, "figures", "output")

GENES = ["aqp4", "aqp1", "aqp2", "aqp5", "aqp6", "aqp8"]
PARALOGUES = GENES[1:]
REF_LEN = dict(aqp1=269, aqp2=271, aqp4=323, aqp5=265, aqp6=282, aqp8=261)
REF = "Homo_sapiens"
QCUT, MIN_ALPHA, MAX_ALPHA, OCC_MIN = 0.20, 0.01, 10.0, 0.50

# ---------------------------------------------------------------- palette --
# Reference categorical palette (light mode), used in fixed slot order.
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
MAGENTA, GREEN, VIOLET, RED = "#e87ba4", "#008300", "#4a3aa7", "#e34948"
SLOTS = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8c8b86", "#e4e3df"
NEUTRAL = "#b9b8b2"                     # "other / not called" marks
SEQ_BLUE = ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
            "#256abf", "#184f95", "#0d366b"]
SEQ_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("seqblue", SEQ_BLUE)

# Residue colours for logos/alignments (grouped by chemistry, 5 groups).
AA_GROUPS = {"DE": RED, "KRH": BLUE, "AVLIMFWC": INK2, "STNQY": GREEN, "PG": YELLOW}
AA_COLOR = {a: c for k, c in AA_GROUPS.items() for a in k}

GENE_LABEL = {g: g.upper() for g in GENES}

# ----------------------------------------------------------------- clades --
CLADE_ORDER = ["Cyclostomata", "Chondrichthyes", "Teleostei: Ostariophysi",
               "Teleostei: Acanthomorpha", "Teleostei: other",
               "Actinopterygii: non-teleost", "Sarcopterygian fishes", "Amphibia",
               "Aves", "Non-avian reptiles", "Monotremata", "Marsupialia", "Eutheria"]
# Display names, harmonised with the manuscript text.
CLADE_LABEL = {
    "Cyclostomata": "Cyclostomata", "Chondrichthyes": "Chondrichthyes",
    "Teleostei: Ostariophysi": "Ostariophysi",
    "Teleostei: Acanthomorpha": "Acanthomorphata",
    "Teleostei: other": "Other teleosts",
    "Actinopterygii: non-teleost": "Non-teleost actinopterygians",
    "Sarcopterygian fishes": "Sarcopterygian fishes", "Amphibia": "Amphibia",
    "Aves": "Aves", "Non-avian reptiles": "Non-avian sauropsids",
    "Monotremata": "Monotremata", "Marsupialia": "Marsupialia", "Eutheria": "Eutheria",
}
# Seven colour groups (the palette validates up to eight slots on adjacent marks).
GROUP_OF = {
    "Eutheria": "Mammalia", "Marsupialia": "Mammalia", "Monotremata": "Mammalia",
    "Aves": "Aves", "Non-avian reptiles": "Non-avian sauropsids",
    "Amphibia": "Amphibia + sarcopterygian fishes",
    "Sarcopterygian fishes": "Amphibia + sarcopterygian fishes",
    "Teleostei: Acanthomorpha": "Acanthomorphata",
    "Teleostei: Ostariophysi": "Other ray-finned fishes",
    "Teleostei: other": "Other ray-finned fishes",
    "Actinopterygii: non-teleost": "Other ray-finned fishes",
    "Chondrichthyes": "Cartilaginous + jawless fishes",
    "Cyclostomata": "Cartilaginous + jawless fishes",
}
GROUP_ORDER = ["Mammalia", "Aves", "Non-avian sauropsids",
               "Amphibia + sarcopterygian fishes", "Acanthomorphata",
               "Other ray-finned fishes", "Cartilaginous + jawless fishes"]
GROUP_COLOR = dict(zip(GROUP_ORDER, SLOTS))
DEPTH_LABEL = {"gnathostomes": "Gnathostomata", "osteichthyans": "Euteleostomi",
               "tetrapods": "Tetrapoda", "mammals": "Mammalia"}


# ------------------------------------------------------------------ style --
def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
        "ytick.color": INK2, "text.color": INK, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "bold", "axes.titlelocation": "left",
        "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
        "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "lines.linewidth": 1.2, "patch.linewidth": 0,
    })


style()


def save(fig, figure, panel):
    """Save one panel as SVG, PDF and PNG in figures/output/<figure>/."""
    d = os.path.join(OUTDIR, figure)
    os.makedirs(d, exist_ok=True)
    stem = os.path.join(d, f"{figure}{panel}")
    for ext in ("svg", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    fig.savefig(f"{stem}.png", dpi=600)
    plt.close(fig)
    print("saved", os.path.relpath(stem, ROOT) + ".{svg,pdf,png}")


def recessive_grid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------- loaders --
def species(gene):
    """Per-species table: accession, clade, set flags, stop, +4..+7, extension."""
    df = pd.read_csv(os.path.join(TAB, f"{gene}_species.tsv"), sep="\t",
                     keep_default_na=False, dtype={"next4": str, "extension": str})
    for c in ("A1", "A2", "A3", "A4", "selection", "orf_codons", "utr_len"):
        df[c] = df[c].astype(int)
    df["group"] = df["clade"].map(GROUP_OF)
    df["context"] = np.where(df["next4"].str.len() == 4,
                             df["stop"].str.replace("T", "U") + "-" +
                             df["next4"].str.replace("T", "U"), "")
    df = df.set_index("species", drop=False)
    df.index.name = None
    return df


def site_table(name):
    """A per-site table written by final_analysis.py, e.g. 'aqp4_joint_mafft_extension'."""
    df = pd.read_csv(os.path.join(TAB, f"{name}.tsv"), sep="\t", keep_default_na=False,
                     na_values=["nan"])
    df["tested"] = df["exclusion_reason"].astype(str).eq("") & df["q"].notna()
    df["purifying"] = df["call"].eq("purifying")
    df["diversifying"] = df["call"].eq("diversifying")
    df["bonf"] = df["bonferroni"].astype(int).eq(1)
    return df


def length_null(gene):
    return pd.read_csv(os.path.join(TAB, f"{gene}_length_null.tsv"), sep="\t")


def length_null_shuffles(gene):
    """Every dinucleotide-shuffle null length behind <gene>_length_null.tsv (200 per species).

    final_analysis.py keeps only each species' null median/mean, so the shuffles are
    regenerated with its own functions and seed (cached in figures/cache/) and checked
    against the table. Returns the table with a `null` column (array of shuffle lengths).
    """
    import random
    import sys
    cache = os.path.join(ROOT, "figures", "cache", f"{gene}_length_null_shuffles.json")
    if not os.path.exists(cache):
        sys.path.insert(0, SCRIPTS)
        cwd = os.getcwd()
        os.chdir(ROOT)
        try:
            import final_analysis as F
            G = F.load_gene(gene)
        finally:
            os.chdir(cwd)
        rec, rng, out = G["rec"], random.Random(F.SEED), {}
        for s in sorted(G["A4"]):          # same order and rng use as F.length_null()
            u = rec[s]["utr"]
            if len(u) < F.NULL_MINLEN:
                continue
            w = u[:F.NULL_WIN]
            if re.search("[^ACGT]", w):
                continue
            cap = len(w) // 3
            out[s] = [F.orf_len(F.dinuc_shuffle(w, rng), 0, cap)[0] for _ in range(F.NULL_REPS)]
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        with open(cache, "w") as fh:
            json.dump(out, fh)
    with open(cache) as fh:
        sh = json.load(fh)
    ln = length_null(gene)
    ln["null"] = [np.array(sh[s]) for s in ln.species]
    med = np.array([np.median(x) for x in ln.null])
    assert len(sh) == len(ln) and np.array_equal(med, ln.null_median.values), \
        f"regenerated shuffles do not reproduce {gene}_length_null.tsv"
    return ln


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def alignment(gene, kind="mafft"):
    """Protein extension alignment (selection set) for a gene/aligner."""
    return read_fasta(os.path.join(RES, gene, f"ext_{kind}_prot.faa"))


ACC = re.compile(r"([A-Za-z]{2}_[\d.]+)")


def raw_utr(gene):
    """3'UTR sequences from the raw NCBI download, keyed by transcript accession."""
    out = {}
    for name, s in read_fasta(os.path.join(DATA, f"3p_utr_{gene}.fna")).items():
        m = ACC.search(name)
        if m:
            out[m.group(1)] = s.upper()
    if not out:   # accession may only be in the description line
        with open(os.path.join(DATA, f"3p_utr_{gene}.fna")) as fh:
            txt = fh.read().split(">")[1:]
        for block in txt:
            head, *seq = block.split("\n")
            m = ACC.search(head)
            if m:
                out[m.group(1)] = "".join(seq).upper()
    return out


def raw_cds(gene):
    out = {}
    with open(os.path.join(DATA, f"cds_{gene}.fna")) as fh:
        txt = fh.read().split(">")[1:]
    for block in txt:
        head, *seq = block.split("\n")
        m = ACC.search(head)
        if m:
            out[m.group(1)] = "".join(seq).upper()
    return out


def fel_json(gene, fit):
    with open(os.path.join(RES, gene, "fel", f"{fit}.json")) as fh:
        return json.load(fh)


def global_omega(gene, fit):
    d = fel_json(gene, fit)
    rd = d["fits"]["Global MG94xREV"]["Rate Distributions"]
    return float(list(rd.values())[0][0][0])


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def fel_sites(gene, fit, nchan=None, occ_min=OCC_MIN):
    """Re-evaluate a HyPhy FEL fit with the manuscript's rules.

    Returns one row per codon column with module ('channel'/'extension'),
    alpha, beta, p, occupancy, tested flag, BH q (within module) and call.
    """
    nchan = REF_LEN[gene] if nchan is None else nchan
    d = fel_json(gene, fit)
    hdr = [h[0] for h in d["MLE"]["headers"]]
    rows = np.array(d["MLE"]["content"]["0"], float)
    aln = read_fasta(os.path.join(RES, gene, "fel", f"{fit}.fna"))
    seqs = list(aln.values())
    ncol = len(seqs[0]) // 3
    occ = np.zeros(ncol)
    for s in seqs:
        for i in range(ncol):
            if s[3 * i:3 * i + 3] != "---":
                occ[i] += 1
    df = pd.DataFrame(dict(column=np.arange(1, ncol + 1),
                           alpha=rows[:, hdr.index("alpha")],
                           beta=rows[:, hdr.index("beta")],
                           p=rows[:, hdr.index("p-value")],
                           occupancy=occ / len(seqs)))
    if nchan == 0:
        df["module"] = "extension"
    else:
        df["module"] = np.where(df["column"] <= nchan, "channel", "extension")
    df["tested"] = ((df["occupancy"] >= occ_min) & (df["alpha"] > MIN_ALPHA)
                    & (df["alpha"] <= MAX_ALPHA))
    df["q"] = np.nan
    for m in ("channel", "extension"):
        sel = df["tested"] & df["module"].eq(m)
        df.loc[sel, "q"] = bh(df.loc[sel, "p"])
    df["purifying"] = df["tested"] & (df["q"] <= QCUT) & (df["beta"] < df["alpha"])
    df["diversifying"] = df["tested"] & (df["q"] <= QCUT) & (df["beta"] > df["alpha"])
    m_ext = (df["tested"] & df["module"].eq("extension")).sum()
    m_ch = (df["tested"] & df["module"].eq("channel")).sum()
    df["bonf"] = df["tested"] & np.where(df["module"].eq("extension"),
                                         df["p"] <= 0.05 / max(m_ext, 1),
                                         df["p"] <= 0.05 / max(m_ch, 1)) & (df["beta"] < df["alpha"])
    df["n_taxa"] = len(seqs)
    return df


def report():
    with open(os.path.join(RES, "final_results.txt")) as fh:
        return fh.read()


def mk_fits():
    """Mk-model fits as written to final_results.txt (section 4)."""
    out = {}
    pat = re.compile(r"^\s+(ER|ARD)\s+lnL (\S+)\s+AIC (\S+)\s+rate gain (\S+) loss (\S+) /My"
                     r"\s+P\(root = present\): flat prior (\S+), stationary prior (\S+)", re.M)
    for m in pat.finditer(report()):
        out[m.group(1)] = dict(lnL=float(m.group(2)), aic=float(m.group(3)),
                               gain=float(m.group(4)), loss=float(m.group(5)),
                               flat=float(m.group(6).rstrip(",")), stat=float(m.group(7)))
    return out


# --------------------------------------------------------------- taxonomy --
_NCBI = None


def lineage(key):
    """NCBI lineage names for a species key (Genus_species), via ete3."""
    global _NCBI
    if _NCBI is None:
        from ete3 import NCBITaxa
        _NCBI = NCBITaxa()
    name = key.replace("_", " ")
    tr = _NCBI.get_name_translator([name, name.split()[0]])
    tid = (tr.get(name) or tr.get(name.split()[0]) or [None])[0]
    if tid is None:
        return []
    lin = _NCBI.get_lineage(tid)
    names = _NCBI.get_taxid_translator(lin)
    return [names[x] for x in lin]


def lca_name(keys):
    lins = [lineage(k) for k in keys]
    lins = [x for x in lins if x]
    if not lins:
        return "?"
    common = lins[0]
    for L in lins[1:]:
        i = 0
        while i < min(len(common), len(L)) and common[i] == L[i]:
            i += 1
        common = common[:i]
    return common[-1] if common else "?"


# ------------------------------------------------------ tree + parsimony --
def cuag_tree():
    """TimeTree species tree for the tree set T, as analysed by asr() in
    final_analysis.py (same function, so the figure cannot diverge from the report).

    Returns (tree, state) with state[tip] = 1 for UGA-CUAG, else 0.
    """
    import sys
    cwd = os.getcwd()
    sys.path.insert(0, SCRIPTS)
    try:
        os.chdir(ROOT)
        import final_analysis as F
        G = F.load_gene("aqp4")
        tax = F.Taxonomy()
        tax.add({s: G["rec"][s]["org"] for s in G["T"]})
        t, state, _ = F.asr_tree(G, tax)
    finally:
        os.chdir(cwd)
    return t, state


def sankoff(t, state):
    cost = {}

    def down(c):
        if c.is_terminal():
            s = state[c.name]
            cost[id(c)] = [0 if s == 0 else math.inf, 0 if s == 1 else math.inf]
        else:
            for ch in c.clades:
                down(ch)
            cost[id(c)] = [sum(min(cost[id(ch)][u] + (u != s) for u in (0, 1))
                               for ch in c.clades) for s in (0, 1)]
    down(t.root)
    return cost


def parsimony_events(t, state, root_state, cost=None, prefer_stay=True):
    """Branch events [(clade_obj, 'gain'|'loss')] for one MP reconstruction."""
    cost = cost or sankoff(t, state)
    events = []

    def up(c, s):
        for ch in c.clades:
            opts = [(cost[id(ch)][u] + (u != s), u) for u in (0, 1)]
            m = min(o[0] for o in opts)
            cands = [u for v, u in opts if v == m]
            u = (s if s in cands else cands[0]) if prefer_stay else \
                ([x for x in cands if x != s] or cands)[0]
            if u != s:
                events.append((ch, "gain" if u == 1 else "loss"))
            up(ch, u)
    up(t.root, root_state)
    return events


def tree_layout(t):
    """x = distance from root (My), y = tip order; returns dicts keyed by id(clade)."""
    x, y = {}, {}
    tips = t.get_terminals()
    for i, tip in enumerate(tips):
        y[id(tip)] = i

    def setx(c, d):
        x[id(c)] = d
        for ch in c.clades:
            setx(ch, d + (ch.branch_length or 0.0))
    setx(t.root, 0.0)

    def sety(c):
        if c.is_terminal():
            return y[id(c)]
        ys = [sety(ch) for ch in c.clades]
        y[id(c)] = (min(ys) + max(ys)) / 2
        return y[id(c)]
    sety(t.root)
    return x, y


def draw_tree(ax, t, state, events=(), ambiguous=(), lw=0.35):
    """Rectangular time tree, tips coloured by CUAG state, events on branches."""
    x, y = tree_layout(t)
    parent = {}
    for c in t.find_clades(order="level"):
        for ch in c.clades:
            parent[id(ch)] = c
    for c in t.find_clades():
        for ch in c.clades:
            ax.plot([x[id(c)], x[id(c)]], [y[id(c)], y[id(ch)]], color=INK2, lw=lw,
                    solid_capstyle="butt")
            ax.plot([x[id(c)], x[id(ch)]], [y[id(ch)], y[id(ch)]], color=INK2, lw=lw,
                    solid_capstyle="butt")
    xmax = max(x.values())
    tips = t.get_terminals()
    ty = [y[id(tp)] for tp in tips]
    tc = [BLUE if state[tp.name] else ORANGE for tp in tips]
    ax.scatter([xmax + 8] * len(tips), ty, c=tc, s=1.2, marker="s", lw=0)

    def mark(ch, kind, hollow=False):
        p = parent[id(ch)]
        xm = (x[id(p)] + x[id(ch)]) / 2
        col = BLUE if kind == "gain" else ORANGE
        ax.scatter([xm], [y[id(ch)]], s=26, marker="^" if kind == "gain" else "v",
                   facecolor="white" if hollow else col, edgecolor=col if hollow else "white",
                   lw=0.9, zorder=5)
        return xm, y[id(ch)]
    pos = []
    for ch, kind in events:
        pos.append((ch, kind, *mark(ch, kind)))
    for ch, kind in ambiguous:
        pos.append((ch, kind, *mark(ch, kind, hollow=True)))
    ax.set_ylim(len(tips) + 2, -3)
    ax.set_xlim(-5, xmax + 14)
    ax.set_yticks([])
    for s in ("left",):
        ax.spines[s].set_visible(False)
    ax.set_xticks(np.arange(0, xmax + 1, 100))
    ax.set_xticklabels([f"{int(xmax - v)}" for v in np.arange(0, xmax + 1, 100)])
    ax.set_xlabel("Millions of years before present")
    return pos, x, y


def clade_strip(ax, t, xpos, width=6):
    """Coloured strip of clade groups next to the tips, with group labels."""
    sp = species("aqp4")
    _, y = tree_layout(t)
    tips = t.get_terminals()
    groups = [GROUP_OF.get(sp.loc[tp.name, "clade"], "") if tp.name in sp.index else ""
              for tp in tips]
    runs, start = [], 0
    for i in range(1, len(groups) + 1):
        if i == len(groups) or groups[i] != groups[start]:
            runs.append((groups[start], start, i - 1))
            start = i
    for g, a, b in runs:
        if not g:
            continue
        ax.add_patch(matplotlib.patches.Rectangle((xpos, a - 0.5), width, b - a + 1,
                                                  color=GROUP_COLOR[g], lw=0))
    # label the largest run of each group
    best = {}
    for g, a, b in runs:
        if g and (g not in best or b - a > best[g][1] - best[g][0]):
            best[g] = (a, b)
    for g, (a, b) in best.items():
        if b - a >= 4:
            ax.text(xpos + width + 3, (a + b) / 2, g, va="center", ha="left", fontsize=6,
                    color=INK2)


def medians_ratio(ln):
    return float(np.median(ln["observed"])) / float(np.median(ln["null_median"]))


def sign_test_one_sided(k, n):
    from scipy.stats import binomtest
    return binomtest(k, n, 0.5, alternative="greater").pvalue
