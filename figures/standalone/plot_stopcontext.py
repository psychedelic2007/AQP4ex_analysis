#!/usr/bin/env python3
"""
Figure plots for the AQP4 stop-context results (Results 1 and 2).

Produces ONE FILE PER PLOT — individual SVG and PDF, no panels — sized for
assembly in Inkscape. Also writes every plotted number to a text file so the
figures can be checked or rebuilt without re-running anything.

Usage
-----
    python3 plot_stopcontext.py --gene aqp4
    python3 plot_stopcontext.py --gene aqp4 --datadir /path/to/fasta
    python3 plot_stopcontext.py --gene aqp4 --only logo,clade

Plots produced
--------------
    logo          sequence logo, stop codon + 4 nt downstream      [MAIN]
    nested        nested conservation: UGA-C -> CUA[G/U] -> CUAG   [MAIN]
    clade         motif frequency by clade, phylogenetic order     [MAIN]
    cladogram     schematic tree with the gain and losses marked   [MAIN]
    contexts      full context frequency ranking                   [SUPP]
    departures    the five species lacking UGA-C                   [SUPP]
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

NT = {"A": "#3aa76d", "C": "#2b6cb0", "G": "#d69e2e", "U": "#c53030", "T": "#c53030"}
C = dict(hi="#1b7f79", lo="#9aa5b1", alt="#8b3a62", warn="#c05746",
         gain="#1b7f79", loss="#c05746", grid="#dfe3e8")


def save(fig, outdir, gene, name):
    for ext in ("svg", "pdf"):
        fig.savefig(os.path.join(outdir, f"{gene}_{name}.{ext}"),
                    bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  [ok] {gene}_{name}.svg / .pdf")


def sanitise(s):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", s)).strip("_")


# ---------------------------------------------------------------- taxonomy --
CLADES = [
    ("Lampreys", {"Petromyzon", "Lampetra", "Lethenteron", "Eptatretus"}),
    ("Cartilaginous fishes", {"Callorhinchus", "Carcharhinus", "Carcharodon",
        "Cetorhinus", "Chiloscyllium", "Amblyraja", "Hemiscyllium", "Hemitrygon",
        "Heptranchias", "Heterodontus", "Hypanus", "Leucoraja", "Mustelus", "Raja",
        "Rhinoraja", "Scyliorhinus", "Stegostoma", "Rhincodon", "Mobula",
        "Pristiophorus", "Torpedo", "Narcine"}),
    ("Basal ray-finned fishes", {"Polypterus", "Erpetoichthys", "Acipenser",
        "Polyodon", "Lepisosteus", "Amia"}),
    ("Other teleosts", {"Anguilla", "Alosa", "Sardina", "Megalops", "Elops",
        "Osmerus", "Hypomesus", "Salmo", "Oncorhynchus", "Esox", "Argentina",
        "Gadus", "Chanos", "Scleropages"}),
    ("Ostariophysi", {"Danio", "Carassius", "Cyprinus", "Abramis", "Barbatula",
        "Garra", "Gobio", "Astyanax", "Colossoma", "Hoplias", "Clarias",
        "Ictalurus", "Silurus", "Pangasianodon", "Pygocentrus", "Salminus",
        "Labeo", "Megalobrama", "Sinocyclocheilus", "Triplophysa", "Misgurnus",
        "Phoxinus", "Leuciscus", "Rutilus", "Scardinius", "Squalius",
        "Onychostoma", "Pimephales", "Puntigrus", "Chanodichthys",
        "Ctenopharyngodon", "Myxocyprinus", "Xyrauchen", "Tachysurus",
        "Hemibagrus", "Paramisgurnus", "Pseudorasbora", "Trichomycterus"}),
    ("Acanthomorpha", {"Takifugu", "Amphiprion", "Acanthopagrus", "Anabas",
        "Ammodytes", "Lampris", "Genypterus", "Oryzias", "Poecilia", "Xiphophorus",
        "Gasterosteus", "Perca", "Sander", "Micropterus", "Seriola", "Scophthalmus",
        "Solea", "Cynoglossus", "Hippocampus", "Syngnathus", "Trachinotus",
        "Boleophthalmus", "Periophthalmus", "Aulonocara", "Rhamphochromis",
        "Copadichromis", "Cyprinodon", "Fundulus", "Kryptolebias", "Nothobranchius",
        "Austrofundulus", "Gambusia", "Poeciliopsis", "Odontesthes", "Leuresthes",
        "Cololabis", "Monopterus", "Mastacembelus", "Channa", "Siniperca", "Lates",
        "Sparus", "Dicentrarchus", "Larimichthys", "Scomber", "Limanda",
        "Platichthys", "Pleuronectes", "Hippoglossus", "Paralichthys",
        "Cyclopterus", "Echiichthys", "Blennius", "Salarias", "Lipophrys",
        "Symphodus", "Labrus", "Myripristis", "Poromitra", "Melanostigma",
        "Gymnodraco", "Trematomus", "Anoplopoma", "Sebastes", "Taurulus",
        "Clinocottus", "Brachionichthys", "Antennarius", "Synchiropus",
        "Festucalex", "Syngnathoides", "Stigmatopora", "Corythoichthys",
        "Phycodurus", "Phyllopteryx", "Vanacampus", "Gouania", "Eucyclogobius",
        "Gobius", "Embiotoca", "Brachistius", "Stegastes", "Acanthochromis",
        "Sphaeramia", "Echeneis", "Betta", "Anarhichas", "Anarrhichthys",
        "Pseudoliparis", "Microchirus", "Trachurus"}),
    ("Lobe-finned fishes", {"Latimeria", "Protopterus", "Neoceratodus"}),
    ("Amphibians", {"Xenopus", "Rana", "Bufo", "Nanorana", "Rhinatrema",
        "Microcaecilia", "Geotrypetes", "Leptobrachium", "Pyxicephalus", "Hyla",
        "Spea", "Ascaphus", "Bombina", "Pelobates", "Engystomops",
        "Eleutherodactylus", "Lithobates", "Dendrobates", "Ranitomeya", "Aquarana",
        "Pelophylax", "Anomaloglossus", "Phyllobates", "Dendropsophus",
        "Mixophyes", "Hyperolius", "Crinia", "Pseudophryne", "Lissotriton",
        "Rhinoderma"}),
    ("Reptiles and birds", {"Gallus", "Meleagris", "Coturnix", "Anas",
        "Taeniopygia", "Anolis", "Python", "Pogona", "Crocodylus", "Chelonia",
        "Chrysemys", "Pelodiscus", "Alligator", "Struthio", "Columba", "Falco",
        "Corvus", "Parus", "Serinus", "Ficedula", "Zonotrichia", "Gavialis",
        "Varanus", "Lacerta", "Podarcis", "Zootoca", "Thamnophis", "Naja",
        "Notechis", "Sceloporus", "Gekko", "Terrapene", "Trachemys", "Mauremys",
        "Dromaius", "Apteryx", "Nothoprocta", "Tinamus", "Aquila", "Haliaeetus",
        "Buteo", "Tyto", "Strix", "Melopsittacus", "Nestor", "Cuculus",
        "Ahaetulla", "Pantherophis", "Pituophis", "Protobothrops", "Pseudonaja",
        "Vipera", "Heloderma", "Eublepharis", "Euleptes", "Paroedura",
        "Sphaerodactylus", "Heteronotia", "Tiliqua", "Rhineura", "Hemicordylus",
        "Natrix", "Natator", "Lepidochelys", "Malaclemys", "Ciconia",
        "Acanthisitta", "Agelaius", "Chamaea", "Colius", "Ammospiza",
        "Antrostomus"}),
]
CLADE_ORDER = [c[0] for c in CLADES] + ["Mammals"]


def clade_of(sp):
    g = sp.split("_")[0]
    for name, members in CLADES:
        if g in members:
            return name
    return "Mammals"


# ------------------------------------------------------------------ loader --
def load(gene, datadir):
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
        if len(u) >= 4:
            out[sp] = dict(acc=acc, stop=s[-3:].replace("T", "U"),
                           down=u[:4].replace("T", "U"), clade=clade_of(sp))
    return out


# -------------------------------------------------------------------- logo --
def letter_patch(ch, x, y, w, h, ax):
    """Draw a single nucleotide letter scaled to (w, h) at (x, y)."""
    tp = TextPath((0, 0), ch, size=1, prop=matplotlib.font_manager.FontProperties(
        family="DejaVu Sans", weight="bold"))
    bb = tp.get_extents()
    if bb.width == 0 or bb.height == 0:
        return
    t = (Affine2D()
         .translate(-bb.x0, -bb.y0)
         .scale(w / bb.width, h / bb.height)
         .translate(x, y))
    ax.add_patch(PathPatch(t.transform_path(tp), facecolor=NT.get(ch, "#666"),
                           edgecolor="none"))


def plot_logo(data, outdir, gene, rep):
    """Information-content logo over the stop codon and the next four bases."""
    seqs = [d["stop"] + d["down"] for d in data.values()]
    n, L = len(seqs), 7
    labels = ["stop\n-3", "stop\n-2", "stop\n-1", "+4", "+5", "+6", "+7"]

    fig, ax = plt.subplots(figsize=(4.2, 2.3))
    lines = []
    for i in range(L):
        col = Counter(s[i] for s in seqs)
        tot = sum(col.values())
        freqs = {b: col.get(b, 0) / tot for b in "ACGU"}
        H = -sum(p * np.log2(p) for p in freqs.values() if p > 0)
        R = 2 - H - (3 / (2 * np.log(2) * tot))          # small-sample correction
        R = max(R, 0)
        y = 0
        for b, p in sorted(freqs.items(), key=lambda t: t[1]):
            if p <= 0.002:
                continue
            h = p * R
            letter_patch(b, i + 0.12, y, 0.76, h, ax)
            y += h
        lines.append(f"pos {labels[i].replace(chr(10),' '):<8} bits={R:5.3f}  " +
                     " ".join(f"{b}={freqs[b]*100:5.1f}%" for b in "ACGU"))

    ax.set_xlim(0, L)
    ax.set_ylim(0, 2.05)
    ax.set_xticks(np.arange(L) + 0.5)
    ax.set_xticklabels(labels)
    ax.set_ylabel("information (bits)")
    ax.axvline(3, color="#333", lw=1.0, ls="--")
    save(fig, outdir, gene, "logo")
    rep.append(("logo  [MAIN]", f"n = {n} species", lines))


# ------------------------------------------------------------------ nested --
def plot_nested(data, outdir, gene, rep):
    """Constraint decaying with distance from the stop codon."""
    n = len(data)
    ugac = sum(1 for d in data.values() if d["stop"] == "UGA" and d["down"][0] == "C")
    cua = sum(1 for d in data.values() if d["stop"] == "UGA" and d["down"][:3] == "CUA")
    cuag = sum(1 for d in data.values() if d["stop"] == "UGA" and d["down"] == "CUAG")

    steps = [("UGA", sum(1 for d in data.values() if d["stop"] == "UGA")),
             ("UGA-C", ugac), ("UGA-CUA", cua), ("UGA-CUAG", cuag)]

    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    y = np.arange(len(steps))[::-1]
    pct = [100 * c / n for _, c in steps]
    ax.barh(y, pct, color=[C["hi"]] * 3 + [C["alt"]], height=.62)
    for yy, (lab, c), p in zip(y, steps, pct):
        ax.text(p + 1.2, yy, f"{p:.1f}%   {c}/{n}", va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels([s[0] for s in steps], fontfamily="monospace")
    ax.set_xlabel("species (%)")
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 25, 50, 75, 100])
    save(fig, outdir, gene, "nested")
    rep.append(("nested  [MAIN]", f"n = {n}",
                [f"{lab:<10}{c:>5} / {n} = {100*c/n:6.2f}%" for lab, c in steps]))


# ---------------------------------------------------------------- contexts --
def plot_contexts(data, outdir, gene, rep):
    """Ranked frequency of every observed stop context."""
    n = len(data)
    ctx = Counter(f"{d['stop']}_{d['down']}" for d in data.values())
    top = ctx.most_common(8)

    fig, ax = plt.subplots(figsize=(3.2, 2.6))
    y = np.arange(len(top))[::-1]
    cols = [C["alt"] if k == "UGA_CUAG" else (C["hi"] if k == "UGA_CUAU" else C["lo"])
            for k, _ in top]
    ax.barh(y, [100 * c / n for _, c in top], color=cols, height=.66)
    for yy, (k, c) in zip(y, top):
        ax.text(100 * c / n + 1.2, yy, f"{100*c/n:.1f}", va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels([k for k, _ in top], fontfamily="monospace", fontsize=7)
    ax.set_xlabel("species (%)")
    ax.set_xlim(0, 100)
    save(fig, outdir, gene, "contexts")
    rep.append(("contexts  [SUPPLEMENTARY]", f"n = {n}",
                [f"{k:<12}{c:>5}  {100*c/n:6.2f}%" for k, c in ctx.most_common(15)]))


# ------------------------------------------------------------------- clade --
def plot_clade(data, outdir, gene, rep):
    """UGA-C and UGA-CUAG frequency by clade, in phylogenetic order."""
    by = defaultdict(list)
    for sp, d in data.items():
        by[d["clade"]].append(d)
    order = [c for c in CLADE_ORDER if c in by]

    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    y = np.arange(len(order))[::-1]
    ugac, cuag, ns = [], [], []
    for c in order:
        v = by[c]
        ns.append(len(v))
        ugac.append(100 * sum(1 for d in v if d["stop"] == "UGA" and d["down"][0] == "C") / len(v))
        cuag.append(100 * sum(1 for d in v if d["stop"] == "UGA" and d["down"] == "CUAG") / len(v))
    ax.barh(y + .18, ugac, height=.34, color=C["lo"], label="UGA-C")
    ax.barh(y - .18, cuag, height=.34, color=C["alt"], label="UGA-CUAG")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c}  (n={m})" for c, m in zip(order, ns)])
    ax.set_xlabel("species (%)")
    ax.set_xlim(0, 108)
    ax.legend(frameon=False, loc="lower right")
    save(fig, outdir, gene, "clade")
    rep.append(("clade  [MAIN]", "grouped by genus assignment; verify boundaries",
                [f"{c:<26}n={m:<5} UGA-C={u:6.1f}%  CUAG={g:6.1f}%"
                 for c, m, u, g in zip(order, ns, ugac, cuag)]))


# --------------------------------------------------------------- cladogram --
def plot_cladogram(data, outdir, gene, rep):
    """Schematic tree with the single gain and the losses marked.

    Topology is drawn, not inferred; it summarises the reconstruction rather
    than reproducing it. Loss counts are read from the data.
    """
    by = defaultdict(list)
    for sp, d in data.items():
        by[d["clade"]].append(d)

    tips = ["Lampreys", "Cartilaginous fishes", "Basal ray-finned fishes",
            "Other teleosts", "Ostariophysi", "Acanthomorpha",
            "Lobe-finned fishes", "Amphibians", "Reptiles and birds", "Mammals"]
    tips = [t for t in tips if t in by]
    pct = {t: 100 * sum(1 for d in by[t] if d["stop"] == "UGA" and d["down"] == "CUAG")
              / len(by[t]) for t in tips}

    depth = {"Lampreys": 0.0, "Cartilaginous fishes": 0.18,
             "Basal ray-finned fishes": 0.40, "Other teleosts": 0.56,
             "Ostariophysi": 0.62, "Acanthomorpha": 0.68,
             "Lobe-finned fishes": 0.44, "Amphibians": 0.60,
             "Reptiles and birds": 0.70, "Mammals": 0.70}

    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    ys = {t: len(tips) - 1 - i for i, t in enumerate(tips)}
    for t in tips:
        x0 = depth.get(t, .5)
        ax.plot([x0, 1.0], [ys[t], ys[t]], color="#555", lw=1.1, zorder=1)
        col = C["hi"] if pct[t] > 50 else C["loss"]
        ax.scatter([1.0], [ys[t]], s=34, color=col, zorder=3, clip_on=False)
        ax.text(1.04, ys[t], f"{t}   {pct[t]:.0f}%", va="center", fontsize=7.5)
    ax.plot([0.0, 0.0], [ys[tips[0]], ys[tips[-1]]], color="#555", lw=1.1)
    ax.scatter([0.14], [np.mean([ys[t] for t in tips[1:]])], marker="*",
               s=150, color=C["gain"], zorder=4)
    ax.text(0.14, np.mean([ys[t] for t in tips[1:]]) + .45,
            "gain of UGA-CUAG\nin the jawed-vertebrate ancestor",
            ha="center", fontsize=7, color=C["gain"])
    for t in tips:
        if pct[t] < 50 and t != "Lampreys":
            ax.scatter([depth.get(t, .5) + .03], [ys[t]], marker="x", s=44,
                       color=C["loss"], zorder=4, lw=1.6)
    ax.set_xlim(-0.04, 1.02)
    ax.set_ylim(-0.8, len(tips) - 0.2)
    ax.axis("off")
    save(fig, outdir, gene, "cladogram")
    rep.append(("cladogram  [MAIN]",
                "topology is schematic; percentages are from the data",
                [f"{t:<26}{pct[t]:6.1f}% CUAG   n={len(by[t])}" for t in tips]))


# -------------------------------------------------------------- departures --
def plot_departures(data, outdir, gene, rep):
    """The species lacking UGA-C, and how each fails."""
    dep = [(sp, d) for sp, d in data.items()
           if not (d["stop"] == "UGA" and d["down"][0] == "C")]
    if not dep:
        return
    dep.sort()
    fig, ax = plt.subplots(figsize=(3.6, 1.4 + .26 * len(dep)))
    y = np.arange(len(dep))[::-1]
    for yy, (sp, d) in zip(y, dep):
        bad_stop = d["stop"] != "UGA"
        bad_p4 = d["down"][0] != "C"
        ax.scatter([0], [yy], s=64, marker="s",
                   color=C["warn"] if bad_stop else C["lo"])
        ax.scatter([1], [yy], s=64, marker="s",
                   color=C["warn"] if bad_p4 else C["lo"])
        ax.text(1.35, yy, f"{sp.replace('_',' ')}   {d['stop']}_{d['down']}",
                va="center", fontsize=7, style="italic")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["stop\ncodon", "+4\nbase"])
    ax.set_xlim(-.5, 4.2)
    ax.set_yticks([])
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    save(fig, outdir, gene, "departures")
    rep.append(("departures  [SUPPLEMENTARY]", f"{len(dep)} species lack UGA-C",
                [f"{sp:<32}{d['stop']}_{d['down']}" for sp, d in dep]))


# -------------------------------------------------------------------- main --
ALL = ["logo", "nested", "clade", "cladogram", "contexts", "departures"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", default="aqp4")
    ap.add_argument("--datadir", default=".")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--only", default=None)
    args = ap.parse_args()

    gene = args.gene.lower()
    outdir = args.outdir or f"figures_{gene}_stopcontext"
    os.makedirs(outdir, exist_ok=True)
    want = set(args.only.split(",")) if args.only else set(ALL)

    data = load(gene, args.datadir)
    print(f"=== {gene}: {len(data)} species -> {outdir}/ ===")
    rep = []

    if "logo" in want:        plot_logo(data, outdir, gene, rep)
    if "nested" in want:      plot_nested(data, outdir, gene, rep)
    if "clade" in want:       plot_clade(data, outdir, gene, rep)
    if "cladogram" in want:   plot_cladogram(data, outdir, gene, rep)
    if "contexts" in want:    plot_contexts(data, outdir, gene, rep)
    if "departures" in want:  plot_departures(data, outdir, gene, rep)

    path = os.path.join(outdir, f"{gene}_stopcontext_numbers.txt")
    with open(path, "w") as fh:
        fh.write(f"NUMBERS BEHIND EVERY PLOT — {gene.upper()} stop context\n")
        fh.write("=" * 70 + "\n")
        fh.write(f"denominator: {len(data)} species with a resolvable stop\n")
        fh.write("codon and at least four downstream nucleotides.\n")
        for name, note, lines in rep:
            fh.write("\n" + "-" * 70 + f"\nPLOT: {name}\nNOTE: {note}\n\n")
            for l in lines:
                fh.write(f"  {l}\n")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
