#!/usr/bin/env python3
"""
omega_profile_ci.py - profile-likelihood 95% CIs for the extension-only global dN/dS
(AQP4 primary fit and the four depth strata; Table 2, Fig. S8b).

scripts/hyphy/FEL_profile.bf is HyPhy 2.5.101 FEL.bf, unchanged up to the global
MG94xREV refit, followed by scripts/hyphy/profile_block.bf: omega is fixed at trial
values with all other parameters (branch lengths, nucleotide rates) re-optimised;
the bounds are found by bracketing plus 12 bisection steps where
logL = logL_max - 1.9207 (chi-square, 1 df, 0.95). FEL site tests are not run.

Inputs are the extension-only alignments and trees written by final_analysis.py
(results_final/aqp4/fel/<fit>.fna/.nwk). One HyPhy log per fit is kept in
results_final/aqp4/omega_ci/ and the table is rebuilt from those logs, so an
existing log is not recomputed unless --force is given.

    python3 scripts/omega_profile_ci.py [--hyphy PATH] [--cpu 6] [--force]
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEL = os.path.join(ROOT, "results_final", "aqp4", "fel")
OUT = os.path.join(ROOT, "results_final", "aqp4", "omega_ci")
BF = os.path.join(ROOT, "scripts", "hyphy", "FEL_profile.bf")
FITS = ["ext_mafft", "depth_gnathostomes_ext", "depth_osteichthyans_ext",
        "depth_tetrapods_ext", "depth_mammals_ext"]


def hyphy_modules(hyphy):
    """<prefix>/share/hyphy/TemplateBatchFiles/SelectionAnalyses/modules for a hyphy binary."""
    prefix = os.path.dirname(os.path.dirname(os.path.realpath(hyphy)))
    d = os.path.join(prefix, "share", "hyphy", "TemplateBatchFiles", "SelectionAnalyses", "modules")
    if not os.path.isdir(d):
        sys.exit(f"HyPhy modules directory not found: {d}")
    return d


def run(fit, hyphy, cpu):
    aln, tree = os.path.join(FEL, f"{fit}.fna"), os.path.join(FEL, f"{fit}.nwk")
    for p in (aln, tree):
        if not os.path.exists(p):
            sys.exit(f"missing {p}; run scripts/final_analysis.py first")
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(BF, tmp)
        os.symlink(hyphy_modules(hyphy), os.path.join(tmp, "modules"))
        log = os.path.join(OUT, f"{fit}.log")
        with open(log + ".tmp", "w") as fh:
            subprocess.run([hyphy, f"CPU={cpu}", "FEL_profile.bf", "--alignment", aln,
                            "--tree", tree, "--output", os.path.join(tmp, "x.json")],
                           cwd=tmp, stdout=fh, stderr=subprocess.STDOUT, check=True)
        os.replace(log + ".tmp", log)


def parse(fit):
    txt = open(os.path.join(OUT, f"{fit}.log")).read()
    mle = re.search(r"PROFILE MLE (\S+) (\S+)", txt)
    ci = re.search(r"PROFILE CI (\S+) (\S+) (\S+)", txt)
    if not (mle and ci):
        sys.exit(f"{fit}.log has no PROFILE result; rerun with --force")
    n = sum(1 for line in open(os.path.join(FEL, f"{fit}.fna")) if line.startswith(">"))
    return n, float(mle.group(1)), float(ci.group(1)), float(ci.group(3)), float(mle.group(2))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hyphy", default=os.environ.get("HYPHY") or shutil.which("hyphy"))
    ap.add_argument("--cpu", type=int, default=6)
    ap.add_argument("--force", action="store_true", help="recompute existing logs")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for fit in FITS:
        if a.force or not os.path.exists(os.path.join(OUT, f"{fit}.log")):
            if not a.hyphy:
                sys.exit("hyphy not found: activate the conda environment or pass --hyphy")
            print(f"profiling {fit} ...", flush=True)
            run(fit, a.hyphy, a.cpu)
    with open(os.path.join(OUT, "omega_profile_ci.tsv"), "w") as fh:
        fh.write("fit\tn_taxa\tomega_mle\tci95_lower\tci95_upper\tlogL_max\n")
        for fit in FITS:
            n, w, lo, hi, L = parse(fit)
            fh.write(f"{fit}\t{n}\t{w:.6f}\t{lo:.6f}\t{hi:.6f}\t{L:.4f}\n")
    print(f"wrote {os.path.relpath(os.path.join(OUT, 'omega_profile_ci.tsv'), ROOT)}")


if __name__ == "__main__":
    main()
