#!/usr/bin/env bash
# run_all.sh - reproduce every number, table and scripted figure panel of the paper.
#
#   bash run_all.sh                  # full pipeline (resumes from results_final/)
#   bash run_all.sh --fresh          # delete results_final/ and recompute everything
#   bash run_all.sh --figures-only   # only redraw figures from the shipped results
#   JOBS=4 THREADS=32 bash run_all.sh
#
# Steps
#   1. check tools and Python packages (see environment.yml)
#   2. NCBI taxonomy for ete3 (downloaded once to ~/.etetoolkit/ if absent)
#   3. scripts/final_analysis.py      species sets, alignments, trees, FEL, report
#   4. scripts/omega_profile_ci.py    profile-likelihood CIs of the extension dN/dS
#   5. figures/run_all.py             every scripted panel -> figures/output/, figures/paper/
#   6. figures/standalone/*.py        Figs S2 and S3
#   7. scripts/fel_to_pymol.py,       PyMOL colouring scripts for the AlphaFold 3 models
#      scripts/seqlogo.py,            (Figs 2f, 3b, 5) -> figures/pymol/, and the Fig. 5
#      figures/standalone/fig5_rates.py rate plots
#
# results_final/ ships with the outputs of the published run. Every external step is
# keyed on the md5 of its inputs, so with the shipped results steps 3-4 only rebuild the
# report (minutes); with --fresh they rerun MAFFT, MUSCLE, PRANK, IQ-TREE and HyPhy
# (many hours, dominated by FEL on the joint alignments).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export PYTHONHASHSEED=0          # fixed set iteration order, so tied entries print in a stable order
JOBS="${JOBS:-4}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"
FRESH=0
FIGURES_ONLY=0
for arg in "$@"; do
    case "$arg" in
        --fresh) FRESH=1 ;;
        --figures-only) FIGURES_ONLY=1 ;;
        -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
        *) echo "unknown option: $arg" >&2; exit 2 ;;
    esac
done

step() { printf '\n=== %s ===\n' "$*"; }

# ------------------------------------------------------------------ 1. checks --
step "1/7 checking tools"
missing=0
need_tools=(python3)
[[ $FIGURES_ONLY -eq 0 ]] && need_tools+=(mafft muscle prank iqtree2 hyphy)
for t in "${need_tools[@]}"; do
    if command -v "$t" >/dev/null 2>&1; then
        printf '  %-8s %s\n' "$t" "$(command -v "$t")"
    else
        printf '  %-8s NOT FOUND\n' "$t"; missing=1
    fi
done
python3 - <<'EOF' || missing=1
import importlib, sys
bad = []
for m in ("Bio", "numpy", "scipy", "pandas", "matplotlib", "logomaker", "ete3"):
    try:
        importlib.import_module(m)
    except Exception:
        bad.append(m)
if bad:
    print("  missing Python packages:", ", ".join(bad))
    sys.exit(1)
print("  Python packages OK")
EOF
if [[ $FIGURES_ONLY -eq 0 ]] && command -v hyphy >/dev/null 2>&1; then
    hv="$(hyphy --version 2>&1 | grep -o 'HYPHY [0-9.]*' | head -1 || true)"
    [[ "$hv" == "HYPHY 2.5.101" ]] || echo "  WARNING: $hv found; the published run used HYPHY 2.5.101 (environment.yml)"
fi
if [[ $missing -ne 0 ]]; then
    echo
    echo "Install the environment first:"
    echo "    conda env create -f environment.yml && conda activate aqp4ex"
    exit 1
fi

# --------------------------------------------------------------- 2. taxonomy --
step "2/7 NCBI taxonomy (ete3)"
python3 - <<'EOF'
import os
from ete3 import NCBITaxa
db = os.path.expanduser("~/.etetoolkit/taxa.sqlite")
print("  using existing", db) if os.path.exists(db) else print("  downloading NCBI taxonomy (once) ...")
NCBITaxa()
EOF

if [[ $FIGURES_ONLY -eq 0 ]]; then
    # ------------------------------------------------------------ 3. analysis --
    if [[ $FRESH -eq 1 ]]; then
        step "removing results_final/ (--fresh)"
        rm -rf results_final
    fi
    step "3/7 analysis (scripts/final_analysis.py, jobs=$JOBS threads=$THREADS)"
    python3 scripts/final_analysis.py --jobs "$JOBS" --threads "$THREADS"

    # ------------------------------------------------------ 4. profile CIs --
    step "4/7 profile-likelihood CIs (scripts/omega_profile_ci.py)"
    python3 scripts/omega_profile_ci.py --cpu "$THREADS"
else
    [[ -f results_final/final_results.txt ]] || { echo "results_final/ is missing; run without --figures-only"; exit 1; }
fi

# ------------------------------------------------------------------ 5. figures --
step "5/7 figure panels (figures/run_all.py)"
python3 figures/run_all.py

# --------------------------------------------------------- 6. Figs S2 and S3 --
step "6/7 Figs S2 and S3 (figures/standalone/)"
SA=figures/output/standalone
python3 figures/standalone/plot_stopcontext.py --gene aqp4 --datadir data --outdir "$SA" --only contexts
python3 figures/standalone/plot_extension.py --gene aqp4 --datadir data --outdir "$SA" --only length_ranked
for ext in svg pdf; do
    cp "$SA/aqp4_contexts.$ext" "figures/paper/FigS2.$ext"
    cp "$SA/aqp4_length_ranked.$ext" "figures/paper/FigS3.$ext"
done

# ------------------------------------------------------------ 7. PyMOL scripts --
step "7/7 PyMOL colouring scripts for the AlphaFold models (figures/pymol/)"
# With AlphaFold models in alphafold/<gene>/{channel,tail}.pdb the scripts load them and the
# residue numbering is checked against the human sequence; without them the scripts colour an
# object named <gene> that you load yourself. Tail models are numbered from 1, i.e. human
# residue - (human protein length + 1).
declare -A HUMAN_LEN=([aqp1]=269 [aqp2]=271 [aqp4]=323 [aqp5]=265 [aqp6]=282 [aqp8]=261)
for g in aqp1 aqp2 aqp4 aqp5 aqp6 aqp8; do
    for region in channel tail; do
        shift=0; [[ $region == tail ]] && shift=$(( -(HUMAN_LEN[$g] + 1) ))
        pdb=(); [[ -f alphafold/$g/$region.pdb ]] && pdb=(--pdb "alphafold/$g/$region.pdb")
        python3 scripts/fel_to_pymol.py --gene "$g" --region "$region" --quantity all \
            --resi-shift "$shift" "${pdb[@]}" --outdir figures/pymol
    done
done
python3 figures/standalone/fig5_rates.py
# conservation of the AQP4 extension for Fig. 2f; offset 0 numbers the 28 residues 1-28, as
# in alphafold/aqp4/tail.pdb
python3 scripts/seqlogo.py --aln results_final/aqp4/ext_mafft_prot.faa --mode ref \
    --reference Homo_sapiens --offset 0 --scale prob --pymol --pymol-object aqp4 \
    --out aqp4_extension_conservation --outdir figures/pymol

step "done"
echo "  report        results_final/final_results.txt"
echo "  tables        results_final/tables/"
echo "  panels        figures/paper/   (manuscript labels; all panels in figures/output/)"
echo "  PyMOL         figures/pymol/   (load the models from alphafold/, see alphafold/README.md)"
