# scripts/

Analysis code. All scripts resolve paths from the repository root, so they can be run from
any directory; `run_all.sh` calls them in the right order.

| Script | Role | Writes |
|---|---|---|
| `final_analysis.py` | The analysis of record: species sets, alignments, gene trees, FEL fits, ancestral states and the report | `results_final/` |
| `omega_profile_ci.py` | Profile-likelihood 95% CIs of the extension-only dN/dS (primary fit and depth strata) | `results_final/aqp4/omega_ci/` |
| `fel_to_pymol.py` | PyMOL scripts colouring the AlphaFold models by site-wise ω, α or β | `figures/pymol/` |
| `seqlogo.py` | Gap-free sequence logos; here used for the conservation colouring of the AQP4 extension | `figures/pymol/` |
| `hyphy/FEL_profile.bf` | HyPhy 2.5.101 `FEL.bf` up to the global MG94×REV fit, followed by `profile_block.bf` | used by `omega_profile_ci.py` |

## final_analysis.py

```bash
python3 scripts/final_analysis.py                    # everything, resumable
python3 scripts/final_analysis.py --jobs 4 --threads 32
python3 scripts/final_analysis.py --stage sets       # species sets only, no external tools
python3 scripts/final_analysis.py --stage report     # rewrite the report from existing outputs
python3 scripts/final_analysis.py --hyphy /path/to/hyphy
```

Stages: `sets` → `align` (MAFFT L-INS-i, MUSCLE, PRANK; IQ-TREE GTR+F+R4 gene trees) →
`fel` (HyPhy FEL) → `report`. Each external job is skipped when its output exists and an
md5 key of its inputs matches, so an interrupted run resumes where it stopped and a changed
input reruns only what depends on it. `--force` rebuilds the derived alignments.

### Species sets (identical rules for every gene)

| Set | Rule |
|---|---|
| A1 | One CDS per species: length a multiple of 3, starts with ATG, ends in a stop codon, no internal stop, protein within 25% of the human length. RefSeq NM_ records are preferred, then the record whose protein length is closest to the human one. |
| A2 | A1 with a 3′-UTR record |
| A3 | A2 whose downstream reading frame ends at an in-frame stop codon inside the 3′-UTR record |
| A4 | A3 without frameshifted records (AQP4 only: frames longer than 35 codons lacking the terminal [A-Z]ET[DEN][LI][LI] motif). For the paralogues A4 = A3. |
| G  | One species per genus from A4: *Homo sapiens* where present, otherwise the longest protein, ties broken alphabetically |
| S  | G with a downstream reading frame of at least 10 codons; used for alignment, gene trees and selection analysis |
| T  | One species per genus from A2 with a resolvable stop context (at least 4 nt after the stop), same rule as G; AQP4 only, used for the ancestral-state analysis |

The genus representative is chosen before the length filter, so the choice is blind to the
trait being studied.

### Selection analysis

Every site-level result comes from a joint channel + extension FEL fit, so α is measured
against the gene-wide synonymous rate. Separate channel-only and extension-only fits give the
module-specific global dN/dS. A codon is tested if at least 50% of sequences have a non-gap
codon in its column and 0.01 < α ≤ 10. p values are corrected by Benjamini-Hochberg within
each gene and module (q ≤ 0.20); Bonferroni is reported alongside.

Robustness runs: MUSCLE and PRANK extension alignments, duplicate extension sequences
removed, and the depth strata Gnathostomata, Euteleostomi, Tetrapoda and Mammalia.

### Ancestral state of UGA-CUAG

The TimeTree tree for set T (`data/species_tree.nwk`) is matched to the species by name, NCBI
species taxid and `data/timetree_synonyms.tsv`; subspecies tips are collapsed after checking
they form a clade. States are reconstructed by unit-cost (Sankoff) parsimony with the root
fixed as present and as absent, and by two-state Mk models (equal rates and all-rates-different)
under flat and stationary root priors.

### Not reported in the paper

Section 6 of the report compares the length of each downstream reading frame with a
dinucleotide-shuffle null (600-nt window, 200 shuffles). It was computed during the
analysis but is not part of the paper; the tables `results_final/tables/<gene>_length_null.tsv`
belong to it.

## omega_profile_ci.py

```bash
python3 scripts/omega_profile_ci.py [--hyphy PATH] [--cpu 6] [--force]
```

For each extension-only fit, ω is fixed at trial values and all other parameters are
re-optimised; the 95% bounds are where the log-likelihood falls 1.92 below its maximum
(bracketing plus 12 bisection steps). One HyPhy log per fit is kept, and the table is rebuilt
from the logs, so existing fits are not recomputed unless `--force` is given.

## fel_to_pymol.py

```bash
python3 scripts/fel_to_pymol.py --gene aqp4 --region tail --resi-shift -324 --pdb alphafold/aqp4/tail.pdb
python3 scripts/fel_to_pymol.py --gene all --region both --quantity all
```

The site tables are rebuilt from the joint FEL fit of record with the functions of
`final_analysis.py`, so the testing rules and calls are the same as in the report. Only tested
sites are coloured; significant sites (q ≤ 0.20) are drawn as sticks. The colour scales are
fixed (log2 ω clipped at ±4; one shared linear scale for α and β), so panels of different
genes compare directly. With `--pdb`, the residue numbering of the model is checked against
the human sequence. Tail models are numbered from 1, so the shift is
−(human protein length + 1): AQP1 −270, AQP2 −272, AQP4 −324, AQP5 −266, AQP6 −283, AQP8 −262.
