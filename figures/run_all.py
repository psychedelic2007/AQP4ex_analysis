"""Regenerate every scripted figure panel of the manuscript and supplement.

Usage:  python3 figures/run_all.py            # all panels
        python3 figures/run_all.py fig2 figS9 # only panels whose script name starts with these

Each panel is produced by its own script in figures/panels/ and saved as SVG, PDF and
PNG (600 dpi) in figures/output/<Figure>/. The script names follow the working numbering
used while the analysis was developed; PAPER below maps them to the panel labels of the
manuscript, and every mapped panel is copied to figures/paper/<label>.{svg,pdf,png}.
Panels that are drawn by hand or rendered in PyMOL are listed in figures/README.md.
"""
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
PAPER_DIR = os.path.join(HERE, "paper")

# manuscript label -> figures/output/<path> (without extension)
PAPER = {
    "Fig1b": "Fig1/Fig1b",
    "Fig1c": "Fig1/Fig1c",
    "Fig1d": "Fig1/Fig1d",
    "Fig2a": "Fig2/Fig2a",
    "Fig2b-d": "Fig2/Fig2c",          # 28-codon logo, 11-codon logo, terminal peptides
    "Fig2e": "FigS6/FigS6a",          # MAFFT extension alignment, representative species
    "Fig3a": "Fig3/Fig3a",
    "Fig3c": "Fig3/Fig3b",
    "Fig3d": "Fig3/Fig3c",
    "Fig3e": "FigS9/FigS9a",          # depth rows (Gn, Eu, Te, Ma) of the robustness panel
    "Fig4a_4e": "Fig4/Fig4a",         # UGA-C share (a) and commonest hexapeptide (e)
    "Fig4b": "FigS15/FigS15b",
    "Fig4c": "Fig4/Fig4c",
    "Fig4d": "Fig4/Fig4d",
    "FigS1": "FigS2/FigS2",
    "FigS4": "FigS4/FigS4b",
    "FigS5a": "FigS7/FigS7a",
    "FigS5b": "FigS7/FigS7c",
    "FigS5c": "FigS7/FigS7d",
    "FigS6": "FigS9/FigS9a",          # MUSCLE, PRANK and duplicate-removed rows
    "FigS7a": "FigS6/FigS6b",
    "FigS7b": "FigS6/FigS6c",
    "FigS8a": "FigS17/FigS17",
    "FigS8b": "FigS9/FigS9b",
    "FigS9": "FigS11/FigS11b",
    "FigS10": "FigS13/FigS13a",
}

scripts = sorted(glob.glob(os.path.join(HERE, "panels", "fig*.py")))
if len(sys.argv) > 1:
    scripts = [s for s in scripts if any(os.path.basename(s).startswith(p) for p in sys.argv[1:])]
failed = []
for s in scripts:
    r = subprocess.run([sys.executable, s], capture_output=True, text=True)
    out = [l for l in r.stdout.splitlines() if l.startswith("saved")]
    print(f"{os.path.basename(s):<14} {'OK' if r.returncode == 0 else 'FAILED'}  {' '.join(out)}")
    if r.returncode:
        failed.append(s)
        print(r.stderr[-2000:])

os.makedirs(PAPER_DIR, exist_ok=True)
for label, src in PAPER.items():
    for ext in ("svg", "pdf", "png"):
        p = os.path.join(OUT, f"{src}.{ext}")
        if os.path.exists(p):
            shutil.copy(p, os.path.join(PAPER_DIR, f"{label}.{ext}"))
print(f"manuscript-labelled copies in {os.path.relpath(PAPER_DIR, os.path.dirname(HERE))}/")
sys.exit(1 if failed else 0)
