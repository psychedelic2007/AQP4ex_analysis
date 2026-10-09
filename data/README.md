# data/

Inputs of the analysis. `scripts/final_analysis.py` reads only these files (plus the NCBI
taxonomy that ete3 keeps in `~/.etetoolkit/`); everything else is derived.

| File | Content |
|---|---|
| `cds_<gene>.fna` | Coding sequences of the vertebrate orthologues, NCBI Orthologs |
| `3p_utr_<gene>.fna` | 3′-UTR records of the same transcripts, NCBI Orthologs |
| `species_tree.nwk` | Time-calibrated species tree from TimeTree 5 for tree set T |
| `timetree_synonyms.tsv` | TimeTree tip names that NCBI taxonomy cannot match by itself |

Genes: AQP4 and the orthodox paralogues AQP1, AQP2, AQP5, AQP6 and AQP8.

| Gene | CDS records | 3′-UTR records |
|---|---|---|
| AQP4 | 1,751 | 1,621 |
| AQP1 | 623 | 597 |
| AQP2 | 488 | 349 |
| AQP5 | 455 | 365 |
| AQP6 | 399 | 325 |
| AQP8 | 746 | 585 |

## NCBI Orthologs downloads

Downloaded on 27 July 2026 from NCBI Orthologs (CDS and 3′-UTR sequence datasets for
the vertebrate orthologues of each human gene). Every 3′-UTR record starts at the nucleotide after the annotated stop codon of
the transcript with the same accession. Headers keep the NCBI format, for example

```
>NM_001308641.1:126-1031 Aqp4 [organism=Mus musculus] [GeneID=11829] [transcript=3] [region=cds]
```

The accession and organism are parsed from the header. The md5 of each file is printed in
section 0 of `results_final/final_results.txt`, so a re-download can be checked against the
published run.

## TimeTree species tree

`species_tree.nwk` was obtained from TimeTree 5 (https://timetree.org) by uploading the list of the 540 species in `results_final/timetree_query.txt`, which
`final_analysis.py` writes on every run (tree set T: one species per genus from A2 with a
resolvable stop context). TimeTree placed 524 of them. It returned 600 tips because some
species come back as several subspecies; these are collapsed after checking that they form
a clade.

Tips are matched to the submitted species by exact name, then by NCBI species taxid
(TimeTree uses synonyms such as *Aquarana catesbeiana* for *Lithobates catesbeianus*), then
through `timetree_synonyms.tsv` (tab-separated: TimeTree tip name, NCBI organism name).
If you change the species sets, upload the new `timetree_query.txt` to TimeTree and replace
this file.
