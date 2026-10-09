#!/usr/bin/env python3
"""
final_analysis.py - one script, from raw NCBI files to every number in the paper.

Everything is rebuilt from the raw inputs. Nothing is read from results_<gene>/,
so the numbers cannot inherit a stale alignment, tree or species set.

INPUTS (data/, next to this scripts/ folder)
    cds_<gene>.fna, 3p_utr_<gene>.fna     NCBI Orthologs downloads
    species_tree.nwk                      TimeTree species tree (ancestral states),
                                          downloaded for results_final/timetree_query.txt
    timetree_synonyms.tsv                 optional: TimeTree tip name <TAB> NCBI organism,
                                          for tips the NCBI taxonomy cannot match by taxid
    ~/.etetoolkit/taxa.sqlite             local NCBI taxonomy (ete3)

TOOLS
    hyphy (>= 2.5.60; the conda 'hyphy' env is used if present), mafft,
    muscle, prank, iqtree2, python3 with biopython, numpy, scipy, ete3

SPECIES SETS (identical rules for every gene)
    A1  one CDS per species: length a multiple of 3, starts ATG, ends in a stop,
        no internal stop, protein within 25% of the human length; NM_ preferred,
        then closest to the human length
    A2  A1 with a 3'UTR record
    A3  A2 whose downstream ORF ends at an in-frame stop inside the record
    A4  A3 minus frameshifted records (AQP4 only: ORF > 35 codons without the
        terminal [A-Z]ET[DEN][LI][LI]? motif). For paralogues A4 = A3.
    S   selection set: one species per genus drawn from A4 (longest protein,
        Homo sapiens protected), then ORF >= 10 codons.
        Genus choice is made before the length filter, so it is blind to the
        trait being studied.
    T   tree set (AQP4 ancestral states): one species per genus drawn from A2
        with a resolvable stop context, same genus rule as S. Written to
        results_final/timetree_query.txt; species_tree.nwk must be TimeTree's
        tree for exactly that list. Tips are matched to T by name, then by NCBI
        taxid (TimeTree uses synonyms), then via timetree_synonyms.tsv.

SELECTION ANALYSES (FEL, HyPhy)
    Every site-level result comes from a JOINT channel+extension fit, so alpha
    is scaled against the gene-wide synonymous rate. Module-specific fits
    (channel-only, extension-only) are run only to report module global dN/dS.
    Testable columns: occupancy >= 50% of taxa (primary; an unfiltered
    sensitivity is also reported), alpha <= 10, alpha > 0.01.
    Benjamini-Hochberg within each gene x module, q <= 0.20; Bonferroni also.

USAGE (paths are resolved from the repository root, so it can be run from anywhere)
    python3 scripts/final_analysis.py                    # everything, resumable
    python3 scripts/final_analysis.py --jobs 4 --threads 32
    python3 scripts/final_analysis.py --stage sets       # fast, no external tools
    python3 scripts/final_analysis.py --stage report     # rewrite report from outputs

OUTPUT
    results_final/final_results.txt              the report
    results_final/<gene>/...                     alignments, trees, FEL JSONs
    results_final/tables/*.tsv                   per-site and per-species tables
"""

import argparse
import hashlib
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO

import numpy as np
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from scipy import optimize, stats

# ----------------------------------------------------------------- settings --
GENES = ["aqp4", "aqp1", "aqp2", "aqp5", "aqp6", "aqp8"]
REF_LEN = dict(aqp1=269, aqp2=271, aqp4=323, aqp5=265, aqp6=282, aqp8=261)
REF = "Homo_sapiens"
LEN_TOL = 0.25            # A1: protein length within 25% of human
LONG_CUT = 35             # AQP4 frameshift inspection threshold (codons)
TERM_MOTIF = re.compile(r"[A-Z]ET[DEN][LI][LI]?$")
MIN_EXT = 10              # selection set: ORF >= 10 codons
MIN_TAXA = 30             # skip selection analysis below this
QCUT = 0.20
MAX_ALPHA = 10.0
MIN_ALPHA = 0.01
OCC_MIN = 0.50            # primary column-occupancy rule
NULL_WIN = 600            # nt of 3'UTR used for the shuffle / other-frame nulls
NULL_MINLEN = 300
NULL_REPS = 200
COMP_WIN = 300            # nt downstream of the ORF stop, for composition
SEED = 20260926
IQ_MODEL = "GTR+F+R4"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root
DATA = os.path.join(ROOT, "data")
OUT = "results_final"     # relative to ROOT; main() and the importers chdir to ROOT


def data(name):
    return os.path.join(DATA, name)


ACC = re.compile(r"([A-Za-z]{2}_[\d.]+)")
ORG = re.compile(r"\[organism=([^\]]+)\]")

CLADE_RULES = [
    ("Cyclostomata", ("Hyperoartia", "Cyclostomata", "Myxini")),
    ("Chondrichthyes", ("Chondrichthyes",)),
    ("Teleostei: Ostariophysi", ("Ostariophysi",)),
    ("Teleostei: Acanthomorpha", ("Acanthomorphata",)),
    ("Teleostei: other", ("Teleostei",)),
    ("Actinopterygii: non-teleost", ("Actinopterygii",)),
    ("Sarcopterygian fishes", ("Actinistia", "Coelacanthimorpha", "Dipnomorpha", "Dipnoi")),
    ("Amphibia", ("Amphibia",)),
    ("Aves", ("Aves",)),
    ("Non-avian reptiles", ("Sauropsida",)),
    ("Monotremata", ("Prototheria",)),
    ("Marsupialia", ("Metatheria",)),
    ("Eutheria", ("Eutheria",)),
]
CLADE_ORDER = [c for c, _ in CLADE_RULES] + ["unassigned"]
DEPTHS = [("gnathostomes", "Gnathostomata"), ("osteichthyans", "Euteleostomi"),
          ("tetrapods", "Tetrapoda"), ("mammals", "Mammalia")]

HYPHY = None   # resolved in main()
LOG = None


# ------------------------------------------------------------------ helpers --
def say(*a):
    msg = " ".join(str(x) for x in a)
    print(msg, flush=True)
    if LOG:
        LOG.write(msg + "\n")
        LOG.flush()


def sanitise(s):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", s)).strip("_")


def read_fasta(p):
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(p, "fasta")}


def write_fasta(d, p, order=None):
    os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        for k in (order or sorted(d)):
            fh.write(f">{k}\n{d[k]}\n")
    os.replace(tmp, p)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(cmd, logfile=None, cwd=None):
    t = time.time()
    with open(logfile or os.devnull, "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=cwd,
                           shell=isinstance(cmd, str))
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): "
                           f"{cmd if isinstance(cmd, str) else ' '.join(cmd)}"
                           f"  -- see {logfile}")
    return time.time() - t


def pct(a, b, d=2):
    return f"{100 * a / b:.{d}f}%" if b else "n/a"


def translate_codon(c):
    if "-" in c or "N" in c or len(c) != 3:
        return None
    return str(Seq(c).translate())


def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    m = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        m = min(m, pvals[i] * n / (rank + 1))
        q[i] = m
    return q


# ----------------------------------------------------------------- taxonomy --
class Taxonomy:
    def __init__(self):
        from ete3 import NCBITaxa
        self.n = NCBITaxa()
        self.lin = {}          # species key -> list of lineage names
        self.unresolved = []

    def add(self, species_orgs):
        names = sorted(set(species_orgs.values()))
        tr = self.n.get_name_translator(names)
        missing = [x for x in names if x not in tr]
        gtr = {}
        if missing:
            gtr = self.n.get_name_translator(sorted({m.split()[0] for m in missing}))
        for key, org in species_orgs.items():
            if key in self.lin:
                continue
            tid = None
            if org in tr:
                tid = tr[org][0]
            elif org.split()[0] in gtr:
                tid = gtr[org.split()[0]][0]
            if tid is None:
                self.lin[key] = []
                self.unresolved.append(org)
                continue
            lineage = self.n.get_lineage(tid)
            names_ = self.n.get_taxid_translator(lineage)
            self.lin[key] = [names_[x] for x in lineage]

    def clade(self, key):
        L = set(self.lin.get(key, []))
        for name, markers in CLADE_RULES:
            if L & set(markers):
                return name
        return "unassigned"

    def within(self, key, taxon):
        return taxon in self.lin.get(key, [])

    def lca_name(self, keys):
        lins = [self.lin.get(k, []) for k in keys if self.lin.get(k)]
        if not lins:
            return "?"
        common = lins[0]
        for L in lins[1:]:
            i = 0
            while i < min(len(common), len(L)) and common[i] == L[i]:
                i += 1
            common = common[:i]
        return common[-1] if common else "?"


# --------------------------------------------------------------- data sets --
def load_gene(gene):
    """Return dict with per-species records and set memberships."""
    utr = {}
    for r in SeqIO.parse(data(f"3p_utr_{gene}.fna"), "fasta"):
        a = ACC.search(r.id) or ACC.search(r.description)
        if a:
            utr[a.group(1)] = str(r.seq).upper()

    n_rec, fail = 0, Counter()
    valid = {}
    for r in SeqIO.parse(data(f"cds_{gene}.fna"), "fasta"):
        n_rec += 1
        a, o = ACC.search(r.id), ORG.search(r.description)
        if not a or not o:
            fail["no accession/organism"] += 1
            continue
        s = str(r.seq).upper()
        if len(s) % 3:
            fail["length not a multiple of 3"] += 1
            continue
        if not s.startswith("ATG"):
            fail["no initiator ATG"] += 1
            continue
        p = str(Seq(s).translate())
        if not p.endswith("*"):
            fail["no terminal stop"] += 1
            continue
        if "*" in p[:-1]:
            fail["internal stop"] += 1
            continue
        valid[a.group(1)] = (o.group(1), s, p[:-1])

    reflen = REF_LEN[gene]
    best = {}
    for acc, (org, s, p) in valid.items():
        sp = sanitise(org)
        k = (not acc.startswith("NM_"), abs(len(p) - reflen))
        if sp not in best or k < best[sp][0]:
            best[sp] = (k, acc, org, s, p)

    rec, outliers = {}, []
    for sp, (_, acc, org, s, p) in best.items():
        if abs(len(p) - reflen) / reflen > LEN_TOL:
            outliers.append((sp, len(p)))
            continue
        u = utr.get(acc, "")
        tail = u[:len(u) // 3 * 3]
        ext = str(Seq(tail).translate(to_stop=True)) if tail else ""
        used = 3 * len(ext)
        stop_ok = bool(u) and (len(tail) >= used + 3 and
                               str(Seq(tail[used:used + 3]).translate()) == "*")
        rec[sp] = dict(sp=sp, acc=acc, org=org, cds=s, prot=p, stop=s[-3:],
                       utr=u, ext=ext, n=len(ext), ext_nt=tail[:used],
                       stop_ok=stop_ok, curated=acc.startswith("NM_"))

    A1 = set(rec)
    A2 = {s for s in A1 if rec[s]["utr"]}
    A3 = {s for s in A2 if rec[s]["stop_ok"]}
    if gene == "aqp4":
        frameshift = sorted(s for s in A3 if rec[s]["n"] > LONG_CUT
                            and not TERM_MOTIF.search(rec[s]["ext"]))
    else:
        frameshift = []
    A4 = A3 - set(frameshift)

    pruned = one_per_genus(A4, rec)
    S = {s for s in pruned if rec[s]["n"] >= MIN_EXT}
    T = one_per_genus({s for s in A2 if context(rec[s]) is not None}, rec)

    return dict(gene=gene, rec=rec, n_records=n_rec, cds_fail=fail,
                n_species_valid=len(best), outliers=sorted(outliers),
                A1=A1, A2=A2, A3=A3, A4=A4, pruned=pruned, S=S, T=T,
                frameshift=frameshift, truncated=sorted(A2 - A3),
                short=sorted(pruned - S, key=lambda x: rec[x]["n"]),
                utr_records=len(utr))


def one_per_genus(species, rec):
    """Genus representatives: Homo sapiens if present, else the longest protein
    (ties alphabetical)."""
    by_genus = defaultdict(list)
    for s in species:
        by_genus[s.split("_")[0]].append(s)
    return {REF if REF in mem else
            sorted(mem, key=lambda x: (-len(rec[x]["prot"]), x))[0]
            for mem in by_genus.values()}


def context(r):
    """(stop codon, next four nt) or None if the UTR is too short."""
    if len(r["utr"]) < 4:
        return None
    return r["stop"], r["utr"][:4]


# --------------------------------------------------- stage: alignment/trees --
def project_codons(prot_aln, nt):
    out = {}
    for k, a in prot_aln.items():
        s = nt[k]
        pieces, i = [], 0
        for c in a:
            if c in "-.":
                pieces.append("---")
            else:
                pieces.append(s[i:i + 3])
                i += 3
        if i != len(s):
            raise RuntimeError(f"codon projection mismatch for {k}: used {i} of {len(s)} nt")
        out[k] = "".join(pieces)
    L = {len(v) for v in out.values()}
    if len(L) != 1:
        raise RuntimeError(f"ragged codon alignment {sorted(L)}")
    return out


def refmask(aln, ref):
    keep = [i for i in range(len(aln[ref]) // 3) if aln[ref][3 * i:3 * i + 3] != "---"]
    return {k: "".join(v[3 * i:3 * i + 3] for i in keep) for k, v in aln.items()}


def tree_write(tree, path):
    """Newick with branch lengths only - no support values or internal labels."""
    def ser(c):
        if c.is_terminal():
            body = c.name
        else:
            body = "(" + ",".join(ser(x) for x in c.clades) + ")"
        bl = c.branch_length if c.branch_length is not None else 0.0
        return f"{body}:{bl:.10g}"
    root = tree.root
    s = "(" + ",".join(ser(x) for x in root.clades) + ");"
    with open(path, "w") as fh:
        fh.write(s + "\n")


def prune_tree(src, keep, dst):
    t = Phylo.read(src, "newick")
    tips = {x.name for x in t.get_terminals()}
    missing = set(keep) - tips
    if missing:
        raise RuntimeError(f"{len(missing)} taxa absent from tree {src}: "
                           f"{sorted(missing)[:5]}")
    for x in [x for x in t.get_terminals() if x.name not in keep]:
        t.prune(x)
    tree_write(t, dst)
    n = len(Phylo.read(dst, "newick").get_terminals())
    if n != len(keep):
        raise RuntimeError(f"pruned tree has {n} tips, expected {len(keep)}")


def align_protein(tool, fasta, out, threads, tree=None, workdir=None):
    if os.path.exists(out) and os.path.getsize(out):
        return
    log = out + ".log"
    if tool == "mafft":
        with open(out + ".tmp", "w") as fh:
            r = subprocess.run(["mafft", "--maxiterate", "1000", "--globalpair",
                                "--anysymbol", "--thread", str(threads), "--quiet",
                                fasta], stdout=fh, stderr=open(log, "w"))
        if r.returncode:
            raise RuntimeError(f"mafft failed on {fasta}")
        os.replace(out + ".tmp", out)
    elif tool == "muscle":
        run(["muscle", "-in", fasta, "-out", out + ".tmp", "-quiet"], log)
        os.replace(out + ".tmp", out)
    elif tool == "prank":
        prefix = out + ".prank"
        cmd = ["prank", f"-d={fasta}", f"-o={prefix}", "-protein", "-f=fasta",
               "-quiet"]
        if tree:
            cmd.insert(2, f"-t={tree}")
        run(cmd, log)
        os.replace(prefix + ".best.fas", out)
    else:
        raise ValueError(tool)


def stage_align(G, threads, force):
    gene, rec, S = G["gene"], G["rec"], sorted(G["S"])
    d = f"{OUT}/{gene}"
    os.makedirs(d, exist_ok=True)
    if len(S) < MIN_TAXA:
        say(f"[{gene}] selection set {len(S)} < {MIN_TAXA}: skipped")
        return
    if REF not in S:
        say(f"[{gene}] WARNING: {REF} not in selection set; residues numbered by column")

    # channel
    write_fasta({s: rec[s]["prot"] for s in S}, f"{d}/channel_prot.faa")
    align_protein("mafft", f"{d}/channel_prot.faa", f"{d}/channel_prot_aln.faa", threads)
    if force or not os.path.exists(f"{d}/channel_codon.fna"):
        pa = read_fasta(f"{d}/channel_prot_aln.faa")
        ca = project_codons(pa, {s: rec[s]["cds"][:-3] for s in S})
        ref = REF if REF in ca else S[0]
        write_fasta(refmask(ca, ref), f"{d}/channel_codon.fna")

    # tree
    tf = f"{d}/tree/{gene}.treefile"
    if not os.path.exists(tf):
        os.makedirs(f"{d}/tree", exist_ok=True)
        say(f"[{gene}] IQ-TREE on {len(S)} taxa")
        run(["iqtree2", "-s", f"{d}/channel_codon.fna", "-m", IQ_MODEL, "-B", "1000",
             "-T", str(threads), "--seed", str(SEED), "--prefix", f"{d}/tree/{gene}",
             "-quiet", "-redo"], f"{d}/tree/iqtree.stdout")
    if not os.path.exists(f"{d}/tree.nwk"):
        prune_tree(tf, set(S), f"{d}/tree.nwk")

    # extension
    write_fasta({s: rec[s]["ext"] for s in S}, f"{d}/ext_prot.faa")
    nt = {s: rec[s]["ext_nt"] for s in S}
    aligners = ["mafft", "muscle", "prank"] if gene == "aqp4" else ["mafft"]
    for tool in aligners:
        pa_path = f"{d}/ext_{tool}_prot.faa"
        align_protein(tool, f"{d}/ext_prot.faa", pa_path, threads,
                      tree=f"{d}/tree.nwk" if tool == "prank" else None)
        cp = f"{d}/ext_{tool}.fna"
        if force or not os.path.exists(cp):
            pa = read_fasta(pa_path)
            pa = {k: v.replace(".", "-") for k, v in pa.items()}
            if set(pa) != set(S):
                raise RuntimeError(f"{tool} alignment taxa differ from selection set")
            write_fasta(project_codons(pa, nt), cp)
    say(f"[{gene}] alignments and tree ready ({len(S)} taxa)")


# ------------------------------------------------------------ FEL job setup --
def build_jobs(G, tax):
    """Write every alignment/tree pair that FEL needs; return job list."""
    gene, rec, S = G["gene"], G["rec"], sorted(G["S"])
    d = f"{OUT}/{gene}"
    if len(S) < MIN_TAXA or not os.path.exists(f"{d}/tree.nwk"):
        return []
    ch = read_fasta(f"{d}/channel_codon.fna")
    jobs = []

    def job(name, aln, taxa, nchan, module):
        taxa = sorted(taxa)
        if len(taxa) < MIN_TAXA:
            say(f"[{gene}] {name}: only {len(taxa)} taxa, skipped")
            return
        jd = f"{d}/fel"
        os.makedirs(jd, exist_ok=True)
        ap, tp = f"{jd}/{name}.fna", f"{jd}/{name}.nwk"
        write_fasta({t: aln[t] for t in taxa}, ap)
        prune_tree(f"{d}/tree.nwk", set(taxa), tp)
        jobs.append(dict(gene=gene, name=name, aln=ap, tree=tp, json=f"{jd}/{name}.json",
                         key=md5(ap) + md5(tp),
                         log=f"{jd}/{name}.log", nchan=nchan, module=module,
                         ntaxa=len(taxa)))

    tools = ["mafft", "muscle", "prank"] if gene == "aqp4" else ["mafft"]
    ext = {t: read_fasta(f"{d}/ext_{t}.fna") for t in tools}
    nchan = len(next(iter(ch.values()))) // 3
    joint = {t: {s: ch[s] + ext[t][s] for s in S} for t in tools}

    for t in tools:
        job(f"joint_{t}", joint[t], S, nchan, "joint")
    job("ext_mafft", ext["mafft"], S, 0, "extension")
    job("channel", ch, S, nchan, "channel")

    # deduplicated (identical extension nucleotide sequences)
    seen = {}
    for s in sorted(S, key=lambda x: (x != REF, x)):
        seen.setdefault(rec[s]["ext_nt"], s)
    dedup = set(seen.values())
    G["dedup_n"] = len(dedup)
    if gene == "aqp4":
        job("joint_dedup", joint["mafft"], dedup, nchan, "joint")
        for label, taxon in DEPTHS:
            sub = {s for s in S if tax.within(s, taxon)}
            job(f"depth_{label}_joint", joint["mafft"], sub, nchan, "joint")
            job(f"depth_{label}_ext", ext["mafft"], sub, 0, "extension")
    else:
        ugac = {s for s in S if (c := context(rec[s])) and c[0] == "TGA"
                and c[1][0] == "C"}
        G["ugac_in_S"] = len(ugac)
        if len(ugac) >= 50:
            job("joint_ugac", joint["mafft"], ugac, nchan, "joint")
    return jobs


def fel_current(j):
    k = j["json"] + ".inputs"
    return (os.path.exists(j["json"]) and os.path.getsize(j["json"]) and
            os.path.exists(k) and open(k).read().strip() == j["key"])


def run_fel(j, cpu):
    if fel_current(j):
        return j, 0.0, "cached"
    cmd = [HYPHY, f"CPU={cpu}", "fel", "--alignment", j["aln"], "--tree", j["tree"],
           "--output", j["json"] + ".tmp"]
    dt = run(cmd, j["log"])
    os.replace(j["json"] + ".tmp", j["json"])
    with open(j["json"] + ".inputs", "w") as fh:
        fh.write(j["key"] + "\n")
    return j, dt, "ran"


# ------------------------------------------------------------ FEL analysis --
def fel_load(path):
    j = json.load(open(path))
    H = [h[0] for h in j["MLE"]["headers"]]
    content = j["MLE"]["content"]
    rows = []
    for k in sorted(content, key=lambda x: int(x)):
        rows += content[k]
    ia, ib, ip = H.index("alpha"), H.index("beta"), H.index("p-value")
    rd = j["fits"]["Global MG94xREV"]["Rate Distributions"]
    omega = None
    for k, v in rd.items():
        if "non-synonymous/synonymous" in k:
            omega = v[0][0]
    return [(r[ia], r[ib], r[ip]) for r in rows], omega, j["input"]["number of sequences"]


def column_info(aln, offset, ref):
    """Per codon column: occupancy, human label, dominant residue."""
    taxa = list(aln)
    L = len(aln[taxa[0]]) // 3
    info = []
    res = offset
    last = offset
    for c in range(L):
        cod = [aln[t][3 * c:3 * c + 3] for t in taxa]
        present = [x for x in cod if x != "---"]
        aa = Counter(a for a in (translate_codon(x) for x in present) if a)
        h = aln[ref][3 * c:3 * c + 3] if ref in aln else "---"
        if h != "---":
            res += 1
            last = res
            label = str(res)
            haa = translate_codon(h) or "?"
        else:
            label = f"ins>{last}"
            haa = "-"
        tot = sum(aa.values())
        dom = aa.most_common(2)
        info.append(dict(col=c + 1, label=label, human=haa, occ=len(present),
                         ntaxa=len(taxa), occf=len(present) / len(taxa),
                         dom=dom[0][0] if dom else "-",
                         domf=dom[0][1] / tot if tot else 0,
                         dom2=f"{dom[1][0]} {100 * dom[1][1] / tot:.1f}%" if len(dom) > 1 else "",
                         humf=(aa.get(haa, 0) / tot if tot and haa not in "-?" else 0),
                         topcodon=Counter(present).most_common(1)[0] if present else ("", 0)))
    return info


def evaluate(rows, info, occ_min):
    """Apply filters, BH and Bonferroni. rows/info are for one module."""
    out = []
    for (a, b, p), ci in zip(rows, info):
        reason = ""
        if ci["occf"] < occ_min:
            reason = f"occupancy {ci['occ']}/{ci['ntaxa']}"
        elif a > MAX_ALPHA:
            reason = f"alpha {a:.1f} > {MAX_ALPHA:g}"
        elif a <= MIN_ALPHA and b <= MIN_ALPHA:
            reason = "alpha and beta ~0"
        elif a <= MIN_ALPHA:
            reason = "alpha at boundary"
        out.append(dict(ci, alpha=a, beta=b, p=p, reason=reason,
                        omega=(b / a if a > 0 else float("nan"))))
    tested = [r for r in out if not r["reason"]]
    for r, q in zip(tested, bh([r["p"] for r in tested])):
        r["q"] = q
        r["bonf"] = r["p"] <= 0.05 / len(tested)
        r["call"] = ("purifying" if q <= QCUT and r["beta"] < r["alpha"] else
                     "diversifying" if q <= QCUT and r["beta"] > r["alpha"] else "-")
    for r in out:
        r.setdefault("q", float("nan"))
        r.setdefault("bonf", False)
        r.setdefault("call", "not tested")
    return out


def analyse_job(j, occ_min=OCC_MIN):
    rows, omega, nseq = fel_load(j["json"])
    aln = read_fasta(j["aln"])
    nchan = j["nchan"]
    res = dict(omega=omega, nseq=nseq, modules={})
    offset = REF_LEN[j["gene"]] + 1
    if j["module"] in ("joint", "channel"):
        ch = {k: v[:3 * nchan] for k, v in aln.items()}
        res["modules"]["channel"] = evaluate(rows[:nchan], column_info(ch, 0, REF), occ_min)
    if j["module"] in ("joint", "extension"):
        ex = {k: v[3 * nchan:] for k, v in aln.items()}
        res["modules"]["extension"] = evaluate(rows[nchan:], column_info(ex, offset, REF),
                                               occ_min)
    return res


def summarise(mod):
    tested = [r for r in mod if r["call"] != "not tested"]
    pur = [r for r in tested if r["call"] == "purifying"]
    div = [r for r in tested if r["call"] == "diversifying"]
    bon = [r for r in tested if r["bonf"] and r["beta"] < r["alpha"]]
    return tested, pur, div, bon


# ------------------------------------------------------- ancestral states --
def species_taxid(n, name):
    """NCBI taxid at species rank for a (sub)species name, trying the binomial
    if the full name is unknown; None if unresolved."""
    for nm in (name, " ".join(name.split()[:2])):
        tid = n.get_name_translator([nm]).get(nm)
        if tid:
            lin = n.get_lineage(tid[0])
            ranks = n.get_rank(lin)
            sp = [x for x in lin if ranks[x] == "species"]
            return sp[-1] if sp else None
    return None


def tree_tip_map(tips, G, tax):
    """Map TimeTree tip names onto species in T: exact name, then the same NCBI
    species (TimeTree returns synonyms, and subspecies of a submitted species),
    then timetree_synonyms.tsv. Returns {tip: species}, {tip: how}."""
    T, rec = G["T"], G["rec"]
    manual = {}
    if os.path.exists(data("timetree_synonyms.tsv")):
        for line in open(data("timetree_synonyms.tsv")):
            if line.strip() and not line.startswith("#"):
                tip, org = line.rstrip("\n").split("\t")[:2]
                manual[sanitise(tip)] = sanitise(org)
    by_sp = defaultdict(list)
    for s in T:
        by_sp[species_taxid(tax.n, rec[s]["org"])].append(s)
    by_sp.pop(None, None)
    m, how = {}, {}
    for t in tips:
        hits = by_sp.get(species_taxid(tax.n, t.replace("_", " ")), [])
        if t in T:
            m[t], how[t] = t, "name"
        elif len(hits) == 1:
            m[t], how[t] = hits[0], "NCBI species"
        elif manual.get(t) in T:
            m[t], how[t] = manual[t], "synonym file"
    return m, how


def asr_tree(G, tax):
    """TimeTree tree with tips renamed to T species and subspecies collapsed,
    and state[species] = 1 for UGA-CUAG. Also used by figures/common.py."""
    rec = G["rec"]
    t = Phylo.read(data("species_tree.nwk"), "newick")
    tips = [x.name for x in t.get_terminals()]
    m, how = tree_tip_map(tips, G, tax)
    extra = sorted(x for x in tips if x not in m)
    if extra:
        raise RuntimeError(f"{len(extra)} tree tips are not in the tree set T, so "
                           f"species_tree.nwk was not built from {OUT}/timetree_query.txt "
                           f"(or needs timetree_synonyms.tsv): {extra[:12]}")
    # TimeTree returns some submitted species as several subspecies tips; a
    # monophyletic group collapses to one tip without changing the species tree
    groups = defaultdict(list)
    for x in tips:
        groups[m[x]].append(x)
    collapsed = {}
    for sp, g in groups.items():
        if len(g) == 1:
            continue
        if len(t.common_ancestor(g).get_terminals()) != len(g):
            raise RuntimeError(f"tips mapping to {sp} are not monophyletic: {g}")
        for x in sorted(g)[1:]:
            t.prune(x)
        collapsed[sp] = g
    unresolved = sorted(G["T"] - set(m.values()))
    for x in t.get_terminals():
        x.name = m[x.name]
    state = {s: 1 if context(rec[s]) == ("TGA", "CTAG") else 0 for s in m.values()}
    return t, state, dict(tips=tips, m=m, how=how, collapsed=collapsed, unresolved=unresolved)


def asr(G, tax, W):
    t, state, info = asr_tree(G, tax)
    tips, m, how = info["tips"], info["m"], info["how"]
    collapsed, unresolved = info["collapsed"], info["unresolved"]
    n1 = sum(state.values())
    W.h("4. ANCESTRAL STATE OF UGA-CUAG  (TimeTree species tree)")
    W(f"tree set T (one species per genus from A2): {len(G['T'])} species submitted to TimeTree")
    W(f"tree tips {len(tips)}, all in T; matched by name {sum(v == 'name' for v in how.values())}, "
      f"by NCBI species {sum(v == 'NCBI species' for v in how.values())}, "
      f"by synonym file {sum(v == 'synonym file' for v in how.values())}")
    for x in sorted(k for k, v in how.items() if v != "name" and m[k] not in collapsed):
        W(f"     {x} -> {m[x]}  ({how[x]})")
    W(f"species returned as several subspecies tips, each a clade collapsed to one tip: "
      f"{len(collapsed)} ({sum(len(g) for g in collapsed.values())} tips)")
    for sp, g in sorted(collapsed.items()):
        W(f"     {sp}: {len(g)} tips")
    W(f"T species not resolved by TimeTree {len(unresolved)}: {', '.join(unresolved)}")
    W(f"tips scored {len(state)}: with UGA-CUAG {n1}, without {len(state) - n1}")
    rootsizes = [len(c.get_terminals()) for c in t.root.clades]
    W(f"root children sizes: {rootsizes}  "
      f"(outgroup = {[tax.lca_name([x.name for x in c.get_terminals()]) for c in t.root.clades]})")

    # Sankoff, unit cost
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
    rc = cost[id(t.root)]
    W(f"parsimony cost with root = absent: {rc[0]:g};  root = present: {rc[1]:g}")
    W("ROOT STATE: " + ("EQUIVOCAL (equal cost)" if rc[0] == rc[1] else
                        f"{'present' if rc[1] < rc[0] else 'absent'} is more parsimonious"))

    def trace(root_state, prefer_stay):
        events = []

        def up(c, s):
            for ch in c.clades:
                opts = [(cost[id(ch)][u] + (u != s), u) for u in (0, 1)]
                m = min(o[0] for o in opts)
                cands = [u for v, u in opts if v == m]
                u = (s if s in cands else cands[0]) if prefer_stay else \
                    ([x for x in cands if x != s] or cands)[0]
                if u != s:
                    tp = [x.name for x in ch.get_terminals()]
                    events.append(("gain" if u == 1 else "loss", len(tp),
                                   tax.lca_name(tp), tp[:3]))
                up(ch, u)
        up(t.root, root_state)
        return events

    for rs in (1, 0):
        if rc[rs] > min(rc):
            continue
        for stay, lab in ((True, "DELTRAN-like"), (False, "ACCTRAN-like")):
            ev = trace(rs, stay)
            W(f"\n  root = {'present' if rs else 'absent'}, {lab}: {len(ev)} changes "
              f"({sum(e[0] == 'gain' for e in ev)} gains, "
              f"{sum(e[0] == 'loss' for e in ev)} losses)")
            for e in sorted(ev, key=lambda e: -e[1]):
                W(f"     {e[0]:<5} {e[1]:>4} tips  {e[2]:<28} e.g. {', '.join(e[3])}")

    # Mk2 maximum likelihood, branch lengths in My
    nodes = []

    def post(c):
        for ch in c.clades:
            post(ch)
        nodes.append(c)
    post(t.root)

    def loglik(a, b, want_root=False):
        s = a + b
        logL = {}
        for c in nodes:
            if c.is_terminal():
                v = np.array([1.0, 0.0]) if state[c.name] == 0 else np.array([0.0, 1.0])
                logL[id(c)] = (v, 0.0)
                continue
            acc = np.ones(2)
            sc = 0.0
            for ch in c.clades:
                bl = max(ch.branch_length or 0.0, 1e-8)
                e = math.exp(-s * bl)
                P = np.array([[b / s + a / s * e, a / s * (1 - e)],
                              [b / s * (1 - e), a / s + b / s * e]])
                v, s_ch = logL[id(ch)]
                acc *= P @ v
                sc += s_ch
            m = acc.max()
            if m <= 0:
                return (-np.inf, None) if want_root else -np.inf
            logL[id(c)] = (acc / m, sc + math.log(m))
        v, sc = logL[id(t.root)]
        if want_root:
            return v, sc
        return sc + math.log(0.5 * v[0] + 0.5 * v[1])

    fits = {}
    for model in ("ER", "ARD"):
        best = None
        for start in (-6, -4, -2):
            if model == "ER":
                f = lambda x: -loglik(math.exp(x[0]), math.exp(x[0]))
                r = optimize.minimize(f, [start], method="Nelder-Mead")
            else:
                f = lambda x: -loglik(math.exp(x[0]), math.exp(x[1]))
                r = optimize.minimize(f, [start, start], method="Nelder-Mead")
            if best is None or r.fun < best.fun:
                best = r
        a = math.exp(best.x[0])
        b = math.exp(best.x[1]) if model == "ARD" else a
        v, _ = loglik(a, b, want_root=True)
        flat = v / v.sum()
        stat = np.array([b / (a + b), a / (a + b)]) * v
        stat = stat / stat.sum()
        k = 1 if model == "ER" else 2
        fits[model] = dict(lnL=-best.fun, a=a, b=b, aic=2 * k + 2 * best.fun,
                           p_flat=flat[1], p_stat=stat[1])
    W("\n  Mk model (maximum likelihood, TimeTree branch lengths in My):")
    for m, f in fits.items():
        W(f"    {m:<4} lnL {f['lnL']:.2f}  AIC {f['aic']:.2f}  rate gain {f['a']:.2e} "
          f"loss {f['b']:.2e} /My   P(root = present): flat prior {f['p_flat']:.3f}, "
          f"stationary prior {f['p_stat']:.3f}")
    W("  NOTE: with a single outgroup lineage, the root state is determined by priors "
      "and branch lengths, not by outgroup evidence.")


# ------------------------------------------------------------- length null --
def orf_len(seq, start, cap_codons):
    n = 0
    for i in range(start, start + 3 * cap_codons, 3):
        c = seq[i:i + 3]
        if len(c) < 3:
            return n, True
        if c in ("TAA", "TAG", "TGA"):
            return n, False
        n += 1
    return n, True


def dinuc_shuffle(s, rng):
    """Altschul-Erickson dinucleotide-preserving shuffle."""
    if len(s) < 3:
        return s
    edges = defaultdict(list)
    for a, b in zip(s[:-1], s[1:]):
        edges[a].append(b)
    last = s[-1]
    verts = [v for v in edges if v != last]
    while True:
        le = {v: rng.randrange(len(edges[v])) for v in verts}
        ok = True
        for v in verts:
            seen, x = set(), v
            while x != last:
                if x in seen:
                    ok = False
                    break
                seen.add(x)
                x = edges[x][le[x]]
            if not ok:
                break
        if ok:
            break
    lists = {}
    for v, e in edges.items():
        e = list(e)
        if v != last:
            keep = e.pop(le[v])
            rng.shuffle(e)
            e.append(keep)
        else:
            rng.shuffle(e)
        lists[v] = e
    ptr = Counter()
    out = [s[0]]
    cur = s[0]
    while ptr[cur] < len(lists.get(cur, [])):
        nxt = lists[cur][ptr[cur]]
        ptr[cur] += 1
        out.append(nxt)
        cur = nxt
    r = "".join(out)
    assert len(r) == len(s), "shuffle length mismatch"
    return r


def geom(p):
    mean = (1 - p) / p
    med = math.ceil(math.log(0.5) / math.log(1 - p)) - 1
    return mean, med


def length_null(G, tax, genus_pruned):
    rec = G["rec"]
    base = G["A4"]
    rng = random.Random(SEED)
    res = dict(n_base=len(base))
    obs = [rec[s]["n"] for s in base]
    res["obs_median"] = float(np.median(obs))
    res["obs_mean"] = float(np.mean(obs))

    # (a) the formula as used in the manuscript: GC of the first 100 nt, A=U, G=C
    gcs = [(rec[s]["utr"][:100].count("G") + rec[s]["utr"][:100].count("C")) / 100
           for s in base if len(rec[s]["utr"]) >= 100]
    gc = float(np.mean(gcs))
    a = (1 - gc) / 2
    g = gc / 2
    p_old = a * a * a + a * a * g + a * g * a
    res["old"] = dict(gc=gc, p=p_old, inv_p=1 / p_old, ratio=res["obs_median"] * p_old)

    # (b) four-base composition downstream of the ORF's own stop
    freqs = []
    for s in base:
        r = rec[s]
        w = r["utr"][3 * r["n"] + 3: 3 * r["n"] + 3 + COMP_WIN]
        w = re.sub("[^ACGT]", "", w)
        if len(w) >= COMP_WIN // 2:
            freqs.append([w.count(x) / len(w) for x in "ACGT"])
    fA, fC, fG, fT = np.mean(freqs, axis=0)
    p_new = fT * fA * fA + fT * fA * fG + fT * fG * fA
    gm, gmed = geom(p_new)
    res["new"] = dict(n=len(freqs), A=fA, C=fC, G=fG, T=fT, p=p_new, mean=gm, median=gmed,
                      ratio_median=res["obs_median"] / gmed if gmed else float("nan"),
                      ratio_mean=res["obs_mean"] / gm)

    # (c) empirical nulls within each species' own 3'UTR
    per = []
    for s in sorted(base):
        u = rec[s]["utr"]
        if len(u) < NULL_MINLEN:
            continue
        w = u[:NULL_WIN]
        if re.search("[^ACGT]", w):
            continue
        cap = len(w) // 3
        o, _ = orf_len(w, 0, cap)
        f1, _ = orf_len(w, 1, cap)
        f2, _ = orf_len(w, 2, cap)
        null = [orf_len(dinuc_shuffle(w, rng), 0, cap)[0] for _ in range(NULL_REPS)]
        per.append(dict(sp=s, obs=o, f1=f1, f2=f2, nmed=float(np.median(null)),
                        nmean=float(np.mean(null)),
                        p=(1 + sum(x >= o for x in null)) / (NULL_REPS + 1),
                        clade=tax.clade(s), pruned=s in genus_pruned))
    res["per"] = per
    return res


def report_null(res, W, gene):
    W(f"  {gene.upper()}  A4 n={res['n_base']}  observed ORF median {res['obs_median']:.1f}, "
      f"mean {res['obs_mean']:.2f} codons")
    o = res["old"]
    W(f"    [manuscript formula] GC(first 100 nt) {100 * o['gc']:.1f}%  P(stop) "
      f"{100 * o['p']:.2f}%  1/P {o['inv_p']:.1f}  median/(1/P) {o['ratio']:.2f}")
    n = res["new"]
    W(f"    [4-base, {COMP_WIN} nt downstream of ORF stop, n={n['n']}] A {n['A']:.3f} "
      f"C {n['C']:.3f} G {n['G']:.3f} U {n['T']:.3f}  P(stop) {100 * n['p']:.2f}%")
    W(f"       geometric null: mean {n['mean']:.1f}, median {n['median']} codons;  "
      f"observed/null  median {n['ratio_median']:.2f}, mean {n['ratio_mean']:.2f}")
    per = res["per"]
    if not per:
        W("    empirical nulls: no species with >= 300 nt clean 3'UTR")
        return
    for lab, sub in (("all A4", per), ("genus-pruned A4", [x for x in per if x["pruned"]])):
        if not sub:
            continue
        obs = [x["obs"] for x in sub]
        nm = [x["nmed"] for x in sub]
        f12 = [max(x["f1"], x["f2"]) for x in sub]
        gt = sum(x["obs"] > x["nmed"] for x in sub)
        lt = sum(x["obs"] < x["nmed"] for x in sub)
        sign = stats.binomtest(gt, gt + lt, 0.5, alternative="greater").pvalue if gt + lt else float("nan")
        try:
            wil = stats.wilcoxon([x["obs"] - max(x["f1"], x["f2"]) for x in sub],
                                 alternative="greater").pvalue
        except ValueError:
            wil = float("nan")
        sig = sum(x["p"] <= 0.05 for x in sub)
        W(f"    [{lab}, n={len(sub)}, first {NULL_WIN} nt, {NULL_REPS} dinucleotide shuffles]")
        W(f"       median observed {np.median(obs):.1f}; median of species' shuffle-null "
          f"medians {np.median(nm):.1f}; ratio {np.median(obs) / max(np.median(nm), 1e-9):.2f}")
        W(f"       species with observed > own null median: {gt}/{gt + lt} (ties excluded), "
          f"one-sided sign test p = {sign:.2e}")
        W(f"       species with empirical P(null >= observed) <= 0.05: {sig}/{len(sub)}")
        W(f"       other frames (+1/+2, same sequence): median of max {np.median(f12):.1f}; "
          f"frame 0 > both in {sum(x['obs'] > max(x['f1'], x['f2']) for x in sub)}/{len(sub)}; "
          f"Wilcoxon frame0 > max(+1,+2) p = {wil:.2e}")
    by = defaultdict(list)
    for x in per:
        by[x["clade"]].append(x)
    rows = []
    for c in CLADE_ORDER:
        sub = by.get(c, [])
        if len(sub) >= 5:
            rows.append((c, len(sub), np.median([x["obs"] for x in sub]),
                         np.median([x["nmed"] for x in sub])))
    k = sum(r[2] > r[3] for r in rows)
    kk = sum(r[2] != r[3] for r in rows)
    W(f"    clade-level (clades with >= 5 species; treats clades, not species, as replicates):")
    for c, nn, a, b in rows:
        W(f"       {c:<30} n={nn:>4}  observed median {a:>5.1f}  null median {b:>5.1f}")
    if kk:
        W(f"       clades with observed > null: {k}/{kk}; one-sided sign test p = "
          f"{stats.binomtest(k, kk, 0.5, alternative='greater').pvalue:.3g}")


# ------------------------------------------------------------------ report --
class Writer:
    def __init__(self):
        self.buf = []

    def __call__(self, s=""):
        self.buf.append(s)

    def h(self, s):
        self.buf += ["", "=" * 100, s, "=" * 100]

    def s(self, s):
        self.buf += ["", s, "-" * len(s)]

    def save(self, p):
        with open(p, "w") as fh:
            fh.write("\n".join(self.buf) + "\n")


def tool_version(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        out = (r.stdout + r.stderr).strip().splitlines()
        return next((l for l in out if l.strip()), "?")[:90]
    except Exception:
        return "not found"


def site_table(W, mod, title, show_all=False):
    W(f"  {title}")
    W(f"  {'col':>4} {'human':>9} {'occ':>9} {'dominant':>14} {'2nd':>10} {'human%':>7}"
      f" {'alpha':>7} {'beta':>7} {'omega':>7} {'p':>9} {'q':>9} {'bonf':>5}  call / reason")
    for r in mod:
        if not show_all and r["occ"] == 0:
            continue
        lab = f"{r['human']}{r['label']}" if r["human"] not in "-?" else r["label"]
        tail = r["call"] if r["call"] != "not tested" else f"[{r['reason']}]"
        q = f"{r['q']:.2e}" if r["q"] == r["q"] else "-"
        W(f"  {r['col']:>4} {lab:>9} {r['occ']:>4}/{r['ntaxa']:<4} "
          f"{r['dom']:>3} {100 * r['domf']:5.1f}%     {r['dom2']:>10} "
          f"{100 * r['humf']:6.1f}% {r['alpha']:7.3f} {r['beta']:7.3f} "
          f"{r['omega']:7.3f} {r['p']:9.2e} {q:>9} {'Y' if r['bonf'] else '':>5}  {tail}")


def write_site_tsv(path, mod):
    with open(path, "w") as fh:
        fh.write("column\thuman_residue\thuman_aa\toccupancy\tn_taxa\tdominant_aa\t"
                 "dominant_freq\thuman_aa_freq\talpha\tbeta\tomega\tp\tq\tbonferroni\t"
                 "call\texclusion_reason\n")
        for r in mod:
            fh.write(f"{r['col']}\t{r['label']}\t{r['human']}\t{r['occ']}\t{r['ntaxa']}\t"
                     f"{r['dom']}\t{r['domf']:.4f}\t{r['humf']:.4f}\t{r['alpha']:.5g}\t"
                     f"{r['beta']:.5g}\t{r['omega']:.5g}\t{r['p']:.5g}\t{r['q']:.5g}\t"
                     f"{int(r['bonf'])}\t{r['call']}\t{r['reason']}\n")


def build_report(genes_data, tax, jobs_by_gene, timings):
    W = Writer()
    os.makedirs(f"{OUT}/tables", exist_ok=True)
    W("FINAL ANALYSIS REPORT - AQP4ex readthrough extension")
    W(f"generated {time.strftime('%Y-%m-%d %H:%M:%S')}   script final_analysis.py "
      f"(md5 {md5(__file__)})")

    W.h("0. PROVENANCE")
    for g in genes_data:
        for f in (f"cds_{g}.fna", f"3p_utr_{g}.fna"):
            W(f"  {f:<22} md5 {md5(data(f))}")
    W(f"  species_tree.nwk       md5 {md5(data('species_tree.nwk'))}")
    W(f"  hyphy   {HYPHY}: {tool_version([HYPHY, '--version'])}")
    W(f"  mafft   {tool_version(['mafft', '--version'])}  (L-INS-i: --maxiterate 1000 --globalpair)")
    W(f"  muscle  {tool_version(['muscle', '-version'])}")
    W(f"  prank   {tool_version(['prank', '-version'])}")
    W(f"  iqtree2 {tool_version(['iqtree2', '--version'])}  model {IQ_MODEL}, UFBoot 1000, seed {SEED}")
    W(f"  settings: A1 length tolerance {LEN_TOL:.0%}; frameshift rule ORF>{LONG_CUT} without "
      f"terminal motif (AQP4 only); selection set ORF>={MIN_EXT}; FEL filters alpha<={MAX_ALPHA:g}, "
      f"alpha>{MIN_ALPHA:g}, column occupancy>={OCC_MIN:.0%} (primary); BH q<={QCUT}; "
      f"shuffle seed {SEED}")
    if tax.unresolved:
        W(f"  taxonomy: {len(tax.unresolved)} organism names unresolved in NCBI: "
          f"{tax.unresolved[:10]}")

    # ---------------------------------------------------------------- sets --
    W.h("1. SPECIES SETS")
    for g, G in genes_data.items():
        rec = G["rec"]
        W.s(f"{g.upper()}")
        W(f"  CDS records {G['n_records']}; failures: "
          + ", ".join(f"{k} {v}" for k, v in G["cds_fail"].items()))
        W(f"  3'UTR records {G['utr_records']}")
        W(f"  species with a valid CDS {G['n_species_valid']}; length outliers (>25%) "
          f"removed {len(G['outliers'])}: {G['outliers'][:8]}")
        W(f"  A1 {len(G['A1'])} | A2 {len(G['A2'])} | A3 {len(G['A3'])} | A4 {len(G['A4'])} | "
          f"genus-pruned A4 {len(G['pruned'])} | SELECTION SET {len(G['S'])}")
        W(f"  curated (NM_) records in A1: {sum(rec[s]['curated'] for s in G['A1'])} "
          f"({pct(sum(rec[s]['curated'] for s in G['A1']), len(G['A1']), 1)})")
        W(f"  truncated (ORF ran off the 3'UTR record) {len(G['truncated'])}: "
          + ", ".join(f"{s}({rec[s]['n']} cod,{len(rec[s]['utr'])} nt)" for s in G["truncated"]))
        if g == "aqp4":
            W(f"  frameshifted {len(G['frameshift'])}: "
              + ", ".join(f"{s}({rec[s]['n']})" for s in G["frameshift"]))
        W(f"  removed from selection set by ORF < {MIN_EXT}: {len(G['short'])}: "
          + ", ".join(f"{s}({rec[s]['n']})" for s in G["short"]))
        if g == "aqp4" and "dedup_n" in G:
            W(f"  deduplicated selection set (identical extension nt): {G['dedup_n']} "
              f"({len(G['S']) - G['dedup_n']} removed)")
        cl = Counter(tax.clade(s) for s in G["S"])
        W("  selection set by clade: " + "; ".join(f"{c} {cl[c]}" for c in CLADE_ORDER if cl[c]))
        with open(f"{OUT}/tables/{g}_species.tsv", "w") as fh:
            fh.write("species\taccession\tclade\tA1\tA2\tA3\tA4\tselection\tstop\tnext4\t"
                     "stop_confirmed\torf_codons\tutr_len\textension\n")
            for s in sorted(G["A1"]):
                r = rec[s]
                c = context(r)
                fh.write(f"{s}\t{r['acc']}\t{tax.clade(s)}\t1\t{int(s in G['A2'])}\t"
                         f"{int(s in G['A3'])}\t{int(s in G['A4'])}\t{int(s in G['S'])}\t"
                         f"{r['stop']}\t{c[1] if c else ''}\t{r['stop_ok']}\t{r['n']}\t"
                         f"{len(r['utr'])}\t{r['ext']}\n")

    G = genes_data["aqp4"]
    rec = G["rec"]

    # -------------------------------------------------------- stop context --
    W.h("2. STOP-CODON CONTEXT (set A2)")
    for g, Gx in genes_data.items():
        ctx = {s: context(Gx["rec"][s]) for s in Gx["A2"]}
        res = {s: c for s, c in ctx.items() if c}
        n = len(res)
        uga = sum(c[0] == "TGA" for c in res.values())
        c4 = sum(c[1][0] == "C" for c in res.values())
        ugac = sum(c[0] == "TGA" and c[1][0] == "C" for c in res.values())
        W(f"  {g.upper()}: A2 {len(Gx['A2'])}, resolvable {n};  UGA {uga} ({pct(uga, n)}), "
          f"+4 C {c4} ({pct(c4, n)}), UGA-C {ugac} ({pct(ugac, n)})")
    ctx = {s: context(rec[s]) for s in G["A2"]}
    ctx = {s: c for s, c in ctx.items() if c}
    n = len(ctx)
    W.s("AQP4 detail")
    cua = sum(c[0] == "TGA" and c[1][:3] == "CTA" and c[1][3] in "GT" for c in ctx.values())
    cuag = sum(c == ("TGA", "CTAG") for c in ctx.values())
    cuau = sum(c == ("TGA", "CTAT") for c in ctx.values())
    W(f"  UGA-CUA[G/U] {cua}/{n} ({pct(cua, n)});  UGA-CUAG {cuag} ({pct(cuag, n)});  "
      f"UGA-CUAU {cuau} ({pct(cuau, n)})")
    W("  all contexts (stop + 4 nt):")
    for (st, nx), k in Counter(ctx.values()).most_common():
        W(f"     {st}-{nx}  {k:>4}  {pct(k, n)}")
    W("  species lacking UGA-C:")
    for s, c in sorted(ctx.items()):
        if not (c[0] == "TGA" and c[1][0] == "C"):
            W(f"     {s:<34} {c[0]}-{c[1]}  {tax.clade(s)}  {rec[s]['acc']}")
    W("  cyclostomes:")
    for s, c in sorted(ctx.items()):
        if tax.clade(s) == "Cyclostomata":
            W(f"     {s:<34} {c[0]}-{c[1]}  {rec[s]['acc']}  curated={rec[s]['curated']}")

    # ------------------------------------------------------------- clades --
    W.h("3. UGA-CUAG BY CLADE (set A2, NCBI taxonomy)")
    W(f"  {'clade':<30}{'n':>5}{'UGA-C':>14}{'CUA[G/U]':>14}{'CUAG':>14}{'CUAU':>14}")
    by = defaultdict(list)
    for s, c in ctx.items():
        by[tax.clade(s)].append(c)
    for cl in CLADE_ORDER:
        L = by.get(cl, [])
        if not L:
            continue
        k = len(L)
        a = sum(c[0] == "TGA" and c[1][0] == "C" for c in L)
        b = sum(c[0] == "TGA" and c[1][:3] == "CTA" and c[1][3] in "GT" for c in L)
        d = sum(c == ("TGA", "CTAG") for c in L)
        e = sum(c == ("TGA", "CTAT") for c in L)
        W(f"  {cl:<30}{k:>5}{a:>6} {pct(a, k, 1):>7}{b:>6} {pct(b, k, 1):>7}"
          f"{d:>6} {pct(d, k, 1):>7}{e:>6} {pct(e, k, 1):>7}")
    unas = sorted(s for s in ctx if tax.clade(s) == "unassigned")
    if unas:
        W(f"  unassigned species: {unas}")

    # ---------------------------------------------------------------- ASR --
    try:
        asr(G, tax, W)
    except Exception as e:
        W(f"  ASR FAILED: {e!r}")
        say(f"WARNING: ancestral-state analysis failed, see section 4 of the report: {e}")

    # ------------------------------------------------ lengths and termini --
    W.h("5. EXTENSION LENGTH AND TERMINAL SEQUENCE (AQP4, set A4)")
    A4 = sorted(G["A4"])
    n4 = len(A4)
    L = Counter(rec[s]["n"] for s in A4)
    W(f"  n = {n4}; median {np.median([rec[s]['n'] for s in A4]):.1f}; "
      f"mean {np.mean([rec[s]['n'] for s in A4]):.2f}")
    for k, v in sorted(L.items(), key=lambda kv: -kv[1]):
        W(f"     {k:>3} codons  {v:>4}  {pct(v, n4)}")
    r25 = sum(25 <= rec[s]["n"] <= 35 for s in A4)
    W(f"  25-35 codons {r25} ({pct(r25, n4)});  > 35 codons: "
      f"{[(s, rec[s]['n']) for s in A4 if rec[s]['n'] > 35]}")
    W(f"  28 + 11 codons: {L[28] + L[11]} ({pct(L[28] + L[11], n4)})")
    t6 = Counter(rec[s]["ext"][-6:] for s in A4 if rec[s]["n"] >= 6)
    t3 = Counter(rec[s]["ext"][-3:] for s in A4 if rec[s]["n"] >= 3)
    n5 = Counter(rec[s]["ext"][:5] for s in A4 if rec[s]["n"] >= 5)
    W(f"  C-terminal hexapeptides (denominator A4 = {n4}):")
    for k, v in t6.most_common(8):
        W(f"     {k}  {v:>4}  {pct(v, n4)}")
    kl = t6["KETDLI"] + t6["KETDIL"]
    W(f"  KETDLI + KETDIL = {kl} ({pct(kl, n4)})")
    dll = sum(bool(re.fullmatch(r"D[LI][LI]", rec[s]["ext"][-3:])) for s in A4)
    W(f"  terminal D-[LI]-[LI]: {dll} ({pct(dll, n4)});  top tripeptides "
      + ", ".join(f"{k} {v}" for k, v in t3.most_common(5)))
    W("  N-terminal pentapeptides:")
    for k, v in n5.most_common(8):
        W(f"     {k}  {v:>4}  {pct(v, n4)}")
    W(f"  ratio KETDLL / top N-pentapeptide = {t6['KETDLL'] / max(n5.most_common(1)[0][1], 1):.2f}")
    W("  11-codon forms:")
    for k, v in Counter(rec[s]["ext"] for s in A4 if rec[s]["n"] == 11).most_common(6):
        W(f"     {k}  {v}")
    W(f"  11-codon forms ending in KETDLL: "
      f"{sum(rec[s]['ext'].endswith('KETDLL') for s in A4 if rec[s]['n'] == 11)}; "
      f"clades: {dict(Counter(tax.clade(s) for s in A4 if rec[s]['n'] == 11))}")
    gc100 = [(rec[s]['utr'][:100].count('G') + rec[s]['utr'][:100].count('C')) / 100
             for s in A4 if len(rec[s]['utr']) >= 100]
    W(f"  mean GC of first 100 nt of 3'UTR: {100 * np.mean(gc100):.2f}% (n={len(gc100)})")

    # --------------------------------------------------------- length null --
    W.h("6. LENGTH OF THE DOWNSTREAM ORF AGAINST CHANCE (all genes, set A4)")
    W("  Observed = codons before the first in-frame stop. Geometric null: mean (1-p)/p, "
      "median ceil(ln .5/ln(1-p))-1.")
    W("  Empirical null: same species' 3'UTR, dinucleotide-shuffled; species are "
      "phylogenetically non-independent, so species-level p-values are liberal.")
    for g, Gx in genes_data.items():
        report_null(Gx["null"], W, g)
        with open(f"{OUT}/tables/{g}_length_null.tsv", "w") as fh:
            fh.write("species\tclade\tobserved\tframe1\tframe2\tnull_median\tnull_mean\t"
                     "empirical_p\tgenus_rep\n")
            for x in Gx["null"]["per"]:
                fh.write(f"{x['sp']}\t{x['clade']}\t{x['obs']}\t{x['f1']}\t{x['f2']}\t"
                         f"{x['nmed']}\t{x['nmean']:.2f}\t{x['p']:.4f}\t{int(x['pruned'])}\n")

    # ------------------------------------------------------------- FEL --
    R = {}
    for g, jobs in jobs_by_gene.items():
        for j in jobs:
            if fel_current(j):
                try:
                    R[(g, j["name"])] = (j, analyse_job(j))
                except Exception as e:
                    say(f"analysis failed {g} {j['name']}: {e!r}")

    def get(g, name):
        return R.get((g, name), (None, None))[1]

    W.h("7. AQP4 SELECTION ANALYSIS - PRIMARY (joint channel+extension, MAFFT, selection set)")
    P = get("aqp4", "joint_mafft")
    if P is None:
        W("  NOT AVAILABLE - FEL has not run")
    else:
        ext = P["modules"]["extension"]
        chn = P["modules"]["channel"]
        write_site_tsv(f"{OUT}/tables/aqp4_joint_mafft_extension.tsv", ext)
        write_site_tsv(f"{OUT}/tables/aqp4_joint_mafft_channel.tsv", chn)
        te, pe, de, be = summarise(ext)
        tc, pc, dc, bc = summarise(chn)
        W(f"  taxa {P['nseq']}; joint global dN/dS {P['omega']:.4f}")
        W(f"  EXTENSION: columns {len(ext)}, tested {len(te)}, purifying (q<={QCUT}) {len(pe)}, "
          f"diversifying {len(de)}, Bonferroni purifying {len(be)}")
        W(f"     expected false discoveries among {len(pe)} calls <= {QCUT * len(pe):.1f}")
        W("     purifying: " + ", ".join(f"{r['human']}{r['label']}(dom {r['dom']})" for r in pe))
        W("     Bonferroni: " + ", ".join(r["label"] for r in be))
        excl = Counter(re.sub(r"[\d./]+", "#", r["reason"]) for r in ext if r["reason"])
        W(f"     excluded columns: {dict(excl)}")
        W(f"  CHANNEL (positive control): tested {len(tc)}, purifying {len(pc)} "
          f"({pct(len(pc), len(tc), 1)}), diversifying {len(dc)} "
          + (f"[{', '.join(r['label'] + ' q=' + format(r['q'], '.3f') for r in dc)}]" if dc else "")
          + f", Bonferroni purifying {len(bc)}")
        site_table(W, ext, "per-column table (extension; occ = non-gap taxa):")
        a_e = [r["alpha"] for r in te]
        a_c = [r["alpha"] for r in tc]
        mw = stats.mannwhitneyu(a_e, a_c, alternative="two-sided")
        rng = np.random.default_rng(SEED)
        pool = np.array(a_e + a_c)
        d0 = np.median(a_e) - np.median(a_c)
        perm = [np.median(x[:len(a_e)]) - np.median(x[len(a_e):])
                for x in (rng.permutation(pool) for _ in range(10000))]
        pp = (1 + sum(abs(x) >= abs(d0) for x in perm)) / 10001
        W(f"\n  SYNONYMOUS RATE, extension vs channel (tested sites, same joint fit):")
        W(f"     median alpha extension {np.median(a_e):.3f} (n={len(a_e)}), channel "
          f"{np.median(a_c):.3f} (n={len(a_c)}); ratio {np.median(a_e) / np.median(a_c):.2f}")
        W(f"     Mann-Whitney two-sided p = {mw.pvalue:.3g}; permutation (10,000) on the "
          f"difference in medians p = {pp:.3g}")
        W(f"     extension sites with alpha < 0.25: "
          + ", ".join(f"{r['label']}(a={r['alpha']:.3f},b={r['beta']:.3f})"
                      for r in te if r["alpha"] < 0.25))
        pos = [r for r in te if r["call"] == "purifying"] or te
        W(f"     among purifying calls: alpha median {np.median([r['alpha'] for r in pos]):.3f} "
          f"range {min(r['alpha'] for r in pos):.3f}-{max(r['alpha'] for r in pos):.3f}; "
          f"beta median {np.median([r['beta'] for r in pos]):.3f} range "
          f"{min(r['beta'] for r in pos):.3f}-{max(r['beta'] for r in pos):.3f}; "
          f"omega median {np.median([r['omega'] for r in pos]):.3f}")
        rho = stats.spearmanr([int(r["label"]) for r in te if r["label"].isdigit()],
                              [r["beta"] for r in te if r["label"].isdigit()])
        W(f"     beta vs human residue position (Spearman): rho {rho.statistic:.2f}, "
          f"p {rho.pvalue:.2g}")
        U = analyse_job(R[("aqp4", "joint_mafft")][0], occ_min=0.0)["modules"]["extension"]
        tu, pu, _, bu = summarise(U)
        W(f"\n  SENSITIVITY - no occupancy filter: tested {len(tu)}, purifying {len(pu)} "
          f"({', '.join(r['label'] for r in pu)}), Bonferroni {len(bu)}")
        for name, lab in (("ext_mafft", "extension-only"), ("channel", "channel-only")):
            x = get("aqp4", name)
            if x:
                W(f"  module global dN/dS, {lab} fit: {x['omega']:.4f}")

    W.h("8. AQP4 ROBUSTNESS")
    runs = [("joint_mafft", "MAFFT (primary)"), ("joint_muscle", "MUSCLE"),
            ("joint_prank", "PRANK"), ("joint_dedup", "MAFFT dedup")]
    runs += [(f"depth_{l}_joint", f"{l} (joint)") for l, _ in DEPTHS]
    calls = {}
    W(f"  {'run':<26}{'taxa':>6}{'dN/dS':>9}{'tested':>8}{'purif.':>8}{'Bonf.':>7}{'div.':>6}"
      f"  ext-only dN/dS")
    for name, lab in runs:
        x = get("aqp4", name)
        if not x:
            W(f"  {lab:<26} not available")
            continue
        t, p, dv, b = summarise(x["modules"]["extension"])
        calls[lab] = {r["label"]: ("B" if r["bonf"] else "q") for r in p}
        eo = ""
        if name.startswith("depth"):
            y = get("aqp4", name.replace("_joint", "_ext"))
            eo = f"{y['omega']:.4f}" if y else "n/a"
        elif name == "joint_mafft":
            y = get("aqp4", "ext_mafft")
            eo = f"{y['omega']:.4f}" if y else "n/a"
        W(f"  {lab:<26}{x['nseq']:>6}{x['omega']:>9.4f}{len(t):>8}{len(p):>8}{len(b):>7}"
          f"{len(dv):>6}  {eo}")
        write_site_tsv(f"{OUT}/tables/aqp4_{name}_extension.tsv", x["modules"]["extension"])
    if calls:
        labs = list(calls)
        allsites = sorted({s for c in calls.values() for s in c},
                          key=lambda s: (int(re.sub(r"\D", "", s) or 0), s))
        W("\n  site x run  (q = q<=0.20, B = also Bonferroni, . = not called)")
        W("  " + f"{'site':>10} " + " ".join(f"{l[:9]:>9}" for l in labs))
        for s in allsites:
            W("  " + f"{s:>10} " + " ".join(f"{calls[l].get(s, '.'):>9}" for l in labs))
        align_runs = [l for l in ["MAFFT (primary)", "MUSCLE", "PRANK", "MAFFT dedup"] if l in calls]
        core = [s for s in allsites if all(s in calls[l] for l in align_runs)]
        coreB = [s for s in allsites if all(calls[l].get(s) == "B" for l in align_runs)]
        W(f"\n  called in all {len(align_runs)} alignment/dedup runs: {core}")
        W(f"  Bonferroni in all {len(align_runs)} alignment/dedup runs: {coreB}")
        every = [s for s in allsites if all(s in calls[l] for l in labs)]
        W(f"  called in every run listed above ({len(labs)}): {every}")

    W.h("9. PARALOGUE COMPARISON (identical rules; joint fits; occupancy >= 50% primary)")
    W(f"  {'gene':<6}{'A2':>5}{'UGA-C':>8}{'A4':>5}{'S':>5}{'ORFmed':>7}{'null':>6}{'top C-hexa':>12}"
      f"{'%A4':>7}{'ext w':>8}{'chan w':>8}{'tested':>7}{'purif':>6}{'%':>6}"
      f"{'| nofilt':>9}{'purif':>6}{'chan%':>7}")
    for g, Gx in genes_data.items():
        rc = Gx["rec"]
        ctx2 = [context(rc[s]) for s in Gx["A2"]]
        ctx2 = [c for c in ctx2 if c]
        ug = sum(c[0] == "TGA" and c[1][0] == "C" for c in ctx2)
        A4g = Gx["A4"]
        t6g = Counter(rc[s]["ext"][-6:] for s in A4g if rc[s]["n"] >= 6)
        top = t6g.most_common(1)[0] if t6g else ("-", 0)
        med = np.median([rc[s]["n"] for s in A4g])
        nl = Gx["null"]["per"]
        nullmed = np.median([x["nmed"] for x in nl]) if nl else float("nan")
        x = get(g, "joint_mafft")
        eo = get(g, "ext_mafft")
        co = get(g, "channel")
        if x:
            t, p, _, _ = summarise(x["modules"]["extension"])
            tc, pc, _, _ = summarise(x["modules"]["channel"])
            U = analyse_job(R[(g, "joint_mafft")][0], occ_min=0.0)["modules"]["extension"]
            tu, pu, _, _ = summarise(U)
            fe = (f"{len(t):>7}{len(p):>6}{pct(len(p), len(t), 1):>6}{len(tu):>9}{len(pu):>6}"
                  f"{pct(len(pc), len(tc), 1):>7}")
        else:
            fe = "   FEL not available"
        W(f"  {g:<6}{len(Gx['A2']):>5}{pct(ug, len(ctx2), 1):>8}{len(A4g):>5}{len(Gx['S']):>5}"
          f"{med:>7.1f}{nullmed:>6.1f}{top[0]:>12}{pct(top[1], len(A4g), 1):>7}"
          f"{(eo['omega'] if eo else float('nan')):>8.3f}"
          f"{(co['omega'] if co else float('nan')):>8.3f}{fe}")
    W("  (ORFmed = median ORF codons in A4; null = median of dinucleotide-shuffle null medians; "
      "ext w / chan w = module-specific global dN/dS; tested/purif = extension, joint fit)")
    for g in genes_data:
        if g == "aqp4":
            continue
        x = get(g, "joint_mafft")
        if not x:
            continue
        ext = x["modules"]["extension"]
        write_site_tsv(f"{OUT}/tables/{g}_joint_mafft_extension.tsv", ext)
        t, p, dv, b = summarise(ext)
        occ = [r for r in ext if r["occ"]]
        W(f"\n  {g.upper()}: extension columns {len(ext)} (with any residue {len(occ)}; "
          f">=50% occupancy {sum(r['occf'] >= OCC_MIN for r in ext)}), taxa {x['nseq']}")
        W(f"     purifying {[r['label'] + '/' + r['dom'] for r in p]}; diversifying "
          f"{[r['label'] for r in dv]}; Bonferroni {len(b)}")
        W(f"     sites beta<alpha among tested: {sum(r['beta'] < r['alpha'] for r in t)}/{len(t)}")
        if t:
            a_e = [r["alpha"] for r in t]
            a_c = [r["alpha"] for r in summarise(x["modules"]["channel"])[0]]
            W(f"     median alpha extension {np.median(a_e):.3f} vs channel {np.median(a_c):.3f}")
        y = get(g, "joint_ugac")
        if y:
            t2, p2, _, _ = summarise(y["modules"]["extension"])
            W(f"     UGA-C species only (n={y['nseq']}): tested {len(t2)}, purifying {len(p2)} "
              f"{[r['label'] for r in p2]}")
        elif "ugac_in_S" in genes_data[g]:
            W(f"     UGA-C species in selection set: {genes_data[g]['ugac_in_S']} "
              f"(< 50, no separate fit)")

    W.h("10. RUN LOG")
    for (g, name), (j, _) in sorted(R.items()):
        W(f"  {g:<5} {name:<24} taxa {j['ntaxa']:>4}  {timings.get((g, name), '')}")
    missing = [(g, j["name"]) for g, js in jobs_by_gene.items() for j in js
               if (g, j["name"]) not in R]
    W(f"  missing/failed FEL runs: {missing if missing else 'none'}")
    W.save(f"{OUT}/final_results.txt")
    say(f"\nreport written: {OUT}/final_results.txt")


# -------------------------------------------------------------------- main --
def main():
    global HYPHY, LOG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", default="all", choices=["all", "sets", "align", "fel", "report"])
    ap.add_argument("--genes", default=",".join(GENES))
    ap.add_argument("--jobs", type=int, default=4, help="concurrent external jobs")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 8)
    ap.add_argument("--hyphy", default=None)
    ap.add_argument("--force", action="store_true", help="rebuild derived alignments")
    a = ap.parse_args()

    os.chdir(ROOT)
    HYPHY = a.hyphy or os.environ.get("HYPHY") or shutil.which("hyphy")
    if not HYPHY and a.stage in ("all", "fel"):
        sys.exit("hyphy not found: activate the conda environment or pass --hyphy /path/to/hyphy")
    os.makedirs(OUT, exist_ok=True)
    LOG = open(f"{OUT}/run.log", "a")
    genes = [g.strip() for g in a.genes.split(",") if g.strip()]
    if "aqp4" not in genes:
        sys.exit("aqp4 must be included")
    for g in genes:
        for f in (f"cds_{g}.fna", f"3p_utr_{g}.fna"):
            if not os.path.exists(data(f)):
                sys.exit(f"missing input data/{f}")
    for f in ("species_tree.nwk",):
        if not os.path.exists(data(f)):
            sys.exit(f"missing input data/{f}")
    say(f"=== final_analysis.py  stage={a.stage}  genes={genes}  hyphy={HYPHY} ===")

    t0 = time.time()
    tax = Taxonomy()
    genes_data = {}
    for g in genes:
        G = load_gene(g)
        tax.add({s: G["rec"][s]["org"] for s in G["A1"]})
        genes_data[g] = G
        say(f"[{g}] A1 {len(G['A1'])} A2 {len(G['A2'])} A3 {len(G['A3'])} "
            f"A4 {len(G['A4'])} S {len(G['S'])}")
    G4 = genes_data["aqp4"]
    with open(f"{OUT}/timetree_query.txt", "w") as fh:
        fh.writelines(f"{G4['rec'][x]['org']}\n" for x in sorted(G4["T"]))
    say(f"[aqp4] tree set T {len(G4['T'])} -> {OUT}/timetree_query.txt")
    tt = Phylo.read(data("species_tree.nwk"), "newick")
    tax.add({x.name: x.name.replace("_", " ") for x in tt.get_terminals()
             if x.name not in tax.lin})
    for g, G in genes_data.items():
        G["null"] = length_null(G, tax, G["pruned"])
    say(f"sets, taxonomy and length nulls done ({time.time() - t0:.0f}s)")

    per_job_threads = max(1, a.threads // a.jobs)
    if a.stage in ("all", "align"):
        with ThreadPoolExecutor(a.jobs) as ex:
            futs = {ex.submit(stage_align, genes_data[g], per_job_threads, a.force): g
                    for g in genes}
            for f in as_completed(futs):
                f.result()

    jobs_by_gene = {g: build_jobs(genes_data[g], tax) for g in genes}
    timings = {}
    if a.stage in ("all", "fel"):
        todo = [j for js in jobs_by_gene.values() for j in js]
        say(f"FEL: {len(todo)} runs, {a.jobs} concurrent x CPU={per_job_threads}")
        with ThreadPoolExecutor(a.jobs) as ex:
            futs = [ex.submit(run_fel, j, per_job_threads) for j in todo]
            for f in as_completed(futs):
                try:
                    j, dt, how = f.result()
                    timings[(j["gene"], j["name"])] = f"{how} {dt / 60:.1f} min"
                    say(f"  FEL {j['gene']} {j['name']}: {how} ({dt / 60:.1f} min)")
                except Exception as e:
                    say(f"  FEL FAILED: {e}")

    build_report(genes_data, tax, jobs_by_gene, timings)
    say(f"total {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
