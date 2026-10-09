# AQP4ex readthrough extension: evolutionary constraint across vertebrates

Code and data for:

> Sangeet, S. *A readthrough-permissive stop context and C-terminal constraint distinguish
> AQP4ex from its orthodox aquaporin paralogues.* bioRxiv (2026). [DOI to be added]

AQP4ex is the C-terminally extended isoform of aquaporin-4 produced by readthrough of
the UGA stop codon. This repository rebuilds every number, table and scripted figure
panel of the paper from the raw NCBI Orthologs downloads:

- the stop-codon context (UGA-C, UGA-CUAG) in AQP4 orthologues from 759 vertebrate species;
- its ancestral state on a TimeTree species tree (parsimony and Mk models);
- the length and terminal sequence of the downstream reading frame;
- site-wise selection (HyPhy FEL) on a joint channel + extension alignment, with
  robustness to aligner, sequence redundancy and phylogenetic depth;
- the same analyses for the orthodox paralogues AQP1, AQP2, AQP5, AQP6 and AQP8.

## Quick start

```bash
git clone <this repository> && cd AQP4ex_github
conda env create -f environment.yml
conda activate aqp4ex
bash run_all.sh
```

`run_all.sh` is the only command you need. It checks the tools, then runs the whole
pipeline in order and prints where each output went.

| Command | What it does | Time |
|---|---|---|
| `bash run_all.sh` | Resumes from the shipped `results_final/`: every external step is keyed on the md5 of its inputs, so only the report, tables and figures are rebuilt | minutes |
| `bash run_all.sh --fresh` | Deletes `results_final/` and recomputes everything (MAFFT, MUSCLE, PRANK, IQ-TREE, HyPhy) | several hours on 32 cores |
| `bash run_all.sh --figures-only` | Redraws the figures from the shipped results; needs only Python | minutes |

Set `JOBS` (concurrent external jobs, default 4) and `THREADS` (default: all cores),
for example `JOBS=4 THREADS=32 bash run_all.sh --fresh`.

## Repository layout

```
AQP4ex_github/
├── run_all.sh              one-shot pipeline (start here)
├── environment.yml         conda environment with the versions of the published run
├── data/                   inputs: NCBI Orthologs FASTA, TimeTree species tree
├── scripts/                analysis code
│   ├── final_analysis.py   species sets -> alignments -> trees -> FEL -> report
│   ├── omega_profile_ci.py profile-likelihood CIs of the extension dN/dS
│   ├── fel_to_pymol.py     PyMOL colouring of the AlphaFold models by FEL rates
│   ├── seqlogo.py          conservation colouring of the AQP4 extension model
│   └── hyphy/              HyPhy batch files for the profile likelihood
├── figures/                figure code
│   ├── run_all.py          draws every scripted panel
│   ├── common.py           shared palette and data loaders
│   ├── panels/             one script per panel
│   ├── standalone/         Figs S2, S3 and the Fig. 5 rate plots
│   ├── paper/              panels under their manuscript labels (generated)
│   └── pymol/              PyMOL scripts and per-site tables (generated)
├── alphafold/              AlphaFold 3 models (add them here, see alphafold/README.md)
└── results_final/          outputs of the published run
    ├── final_results.txt   the report: every number in the paper, by section
    ├── tables/             per-species and per-site tables
    └── <gene>/             alignments, gene trees, FEL fits
```

Each folder has its own README with details.

## Where to find a number from the paper

Every value in the manuscript is in `results_final/final_results.txt`:

| Section | Content | Paper |
|---|---|---|
| 0 | Provenance: md5 of every input, tool versions | Methods |
| 1 | Species sets, with the species removed at each step | Table 1 |
| 2, 3 | Stop-codon context, overall and by clade | Fig. 1b, 1c; Figs S1, S2 |
| 4 | Ancestral state of UGA-CUAG (parsimony, Mk models) | Fig. 1d |
| 5 | Reading-frame length and terminal sequence | Fig. 2; Figs S3, S4 |
| 6 | Reading-frame length against a dinucleotide-shuffle null | computed, not reported in the paper |
| 7 | AQP4 site-wise selection, primary analysis | Fig. 3; Fig. S5 |
| 8 | Robustness: aligners, duplicates removed, depth strata | Table 2; Fig. 3e; Figs S6-S8 |
| 9 | Paralogue comparison | Table 3, Table 4; Figs 4, 5; Figs S8-S10 |
| 10 | Run log of every external job | |

The species sets are defined in
`scripts/README.md` and in Table 1 of the paper:

| Set | Definition | AQP4 |
|---|---|---|
| A1 | One CDS per species meeting the sequence criteria | 759 |
| A2 | A1 with a 3′-UTR record | 697 |
| A3 | A2 whose downstream reading frame ends at an in-frame stop inside the record | 687 |
| A4 | A3 without the two frameshifted records | 685 |
| G  | One species per genus from A4 | 530 |
| S  | G with a downstream reading frame of at least 10 codons (selection set) | 524 |
| T  | One species per genus from A2 with a resolvable stop context (tree set) | 540 submitted, 524 placed |

## Notes on reproducibility

- **HyPhy version.** FEL results depend on the HyPhy version. The published run used
  HyPhy 2.5.101; `run_all.sh` warns if a different version is first on your `PATH`.
- **NCBI taxonomy.** Clade assignment and TimeTree name matching use the NCBI taxonomy
  through ete3. The published run used a copy downloaded on 30 July 2026. On first use ete3
  downloads the current taxonomy to `~/.etetoolkit/`; later renamings can move a few
  species between clades or change how TimeTree synonyms are matched.
- **TimeTree tree.** `data/species_tree.nwk` was downloaded from TimeTree 5
  (timetree.org) for the species listed in `results_final/timetree_query.txt`
  (tree set T). `final_analysis.py` rewrites that list on every run; if the tree contains
  species that are not in the list, the ancestral-state analysis stops with an error and the
  report says so. See `data/README.md`.
- **Deterministic output.** `run_all.sh` sets `PYTHONHASHSEED=0` so that tied entries are
  always listed in the same order. Random steps (IQ-TREE, shuffles) use fixed seeds.

## Citation

If you use this code or data, please cite the paper above. [Citation and DOI to be added]

## Licence

[Add a licence before making the repository public, e.g. MIT for the code and CC BY 4.0
for the data and results.]

## Contact

Satyam Sangeet, School of Physics, The University of Sydney (satyam.sangeet@sydney.edu.au)
