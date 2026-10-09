# results_final/

Outputs of the published run of `scripts/final_analysis.py`. They are included so the
report and figures can be rebuilt without rerunning the alignments and FEL fits:
`bash run_all.sh` reuses them, `bash run_all.sh --fresh` deletes and recomputes them.

| Path | Content |
|---|---|
| `final_results.txt` | The report: every number in the paper, by section (see the main README) |
| `timetree_query.txt` | The 540 species of tree set T, as submitted to TimeTree |
| `tables/<gene>_species.tsv` | One row per species: accession, clade, membership of A1-A4 and S, stop codon, the four following nucleotides, reading-frame length, 3′-UTR length and the extension sequence |
| `tables/<gene>_joint_mafft_extension.tsv` | Site-wise FEL results for the downstream frame (α, β, ω, p, q, call) |
| `tables/aqp4_joint_mafft_channel.tsv` | The same for the AQP4 channel domain |
| `tables/aqp4_joint_{muscle,prank,dedup}_extension.tsv` | Robustness runs |
| `tables/aqp4_depth_<stratum>_joint_extension.tsv` | Depth strata |
| `tables/<gene>_length_null.tsv` | Shuffle-null lengths (computed, not reported in the paper) |
| `<gene>/channel_*.f*`, `<gene>/ext_*.f*` | Protein and codon alignments of the selection set |
| `<gene>/tree.nwk`, `<gene>/tree/` | IQ-TREE gene tree (GTR+F+R4, 1,000 ultrafast bootstraps) |
| `<gene>/fel/<fit>.{fna,nwk,json}` | Input alignment, pruned tree and HyPhy FEL output of every fit |
| `aqp4/omega_ci/` | Profile-likelihood CIs of the extension dN/dS and the HyPhy logs behind them |

FEL fits per gene: `channel`, `ext_mafft` (module-specific dN/dS) and `joint_mafft`
(site-level results). AQP4 also has `joint_muscle`, `joint_prank`, `joint_dedup` and
`depth_<stratum>_{joint,ext}`; AQP5 has `joint_ugac` (the 124 species of S that carry UGA-C).
