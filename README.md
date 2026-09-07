# *Fusarium oxysporum* MRL8996 — Genome assembly & structural effectorome

Command-line strings, custom scripts, and pipelines used for the chromosome-level
genome assembly of *Fusarium oxysporum* strain MRL8996 and the structural
characterisation of its predicted effectorome.

## Overview

- **Assembly:** Oxford Nanopore long reads (Flye) + Hi-C scaffolding (Juicer / 3D-DNA)
- **Total length:** `52.32` Mb
- **Chromosome-scale pseudomolecules:** `16`
- **Final effectorome:** `513` proteins (SignalP ∩ EffectorP, transmembrane-filtered)
- **Structural analyses:** all-vs-all TM-align + DALI, Louvain structural families, AMAPEC antimicrobial-activity prediction

## Repository structure

```
MRL8996_Genome_Project/
├── 01_Genome_Assembly/
│   ├── 01_flye_assembly.sh            # ONT de novo assembly (Flye)
│   ├── 02_hic_scaffolding_juicer.sh   # Hi-C scaffolding (Juicer + 3D-DNA)
│   └── 03_edta_te_annotation.sh       # transposable element annotation (EDTA)
├── 02_Effector_Prediction/
│   ├── 01_secretome_pipeline.sh       # SignalP 6 secretome + mature-seq extraction
│   ├── 02_effectorp_run.sh            # EffectorP 3.0 + intersect + DeepTMHMM
│   ├── 03_parse_effectors.py          # drop TM-containing -> final effectorome
│   ├── 04_structure_prediction.sh     # ColabFold/ESMFold (MATURE seqs) -> PDBs
│   └── 05_amapec_prediction.sh        # AMAPEC antimicrobial-activity prediction
├── 03_Structural_Network/
│   ├── 00_dali_allvsall.sh            # DaliLite v5 all-vs-all -> 'ordered' matrix
│   ├── 01_network_construction.R      # matrix -> igraph -> Louvain families -> figure
│   ├── 02_network_visualization.py    # Python render + Cytoscape GraphML export
│   └── 03_family_representatives.R    # per-family medoid + PyMOL render script
└── 04_Heatmap_Analysis/
    ├── 01_expression_matrix_prep.sh   # TM-align all-vs-all similarity matrix
    └── 02_plot_heatmap.R              # ComplexHeatmap + AMAPEC annotation
```

## Pipeline notes (important for reproducibility)

- **Mature sequences for structure prediction.** AMAPEC explicitly requires structures
  predicted from mature sequences; folding with the signal peptide distorts geometry
  and biases the antimicrobial-activity prediction.
- **AMAPEC output path.** AMAPEC v1.0 writes `AMAPEC_Effectorome/prediction.csv`,
  but `02_plot_heatmap.R` currently reads `AMAPEC_Effectorome/results/prediction.csv`.
  Reconcile the two before plotting.
- **TM-score length bias.** `Max_TM_Score` is normalised by the shorter chain, which
  inflates similarity for size-mismatched pairs — relevant when interpreting edges
  between very different-length effectors.
- **DALI is native all-vs-all** (`dali.pl --matrix`), not snakedali (query-vs-database).

## Software

Flye, BWA, samtools, Juicer, 3D-DNA, Juicebox; EDTA (via Docker); SignalP 6,
EffectorP 3.0, DeepTMHMM (via BioLib); ColabFold / ESMFold; DaliLite v5, TM-align;
AMAPEC; R (igraph, ggraph, ComplexHeatmap, circlize, tidyverse), Python
(networkx, matplotlib, pandas).

## Citation

> Doddi, A.; Puebla Planas, G. *et al.* (2026). A near-complete genome assembly of *Fusarium oxysporum* keratitis  isolate MRL8996.

Tool citations: Flye (Kolmogorov et al. 2019); Juicer (Durand et al. 2016);
3D-DNA (Dudchenko et al. 2017); EDTA (Ou et al. 2019); SignalP 6.0 (Teufel et al. 2022); EffectorP 3.0
(Sperschneider & Dodds 2022); DeepTMHMM (Hallgren et al. 2022); ColabFold
(Mirdita et al. 2022); DALI (Holm 2022); TM-align (Zhang & Skolnick 2005);
AMAPEC (Mesny & Thomma 2024, doi:10.1101/2024.01.04.574150).

## Data availability

Final assembly: NCBI GenBank `GCA_009746015.2`. Raw ONT and Hi-C reads: SRX34109262 & SRX34109263
BioProject `PRJNA554890`. Large sequence files are archived at NCBI, not in this repo.

## 05_Revision_2026 - analyses added during revision

Scripts written to answer the reviewers of the Scientific Data submission. They
run in numerical order and are documented one by one in
`05_Revision_2026/README.md`, which gives for each script its inputs, its
outputs and the figure or table of the manuscript it produces. Paths inside the
scripts are absolute to the analysis machine and have to be adapted.

| Scripts | Analysis |
|---|---|
| 01-04 | Three-way comparison of the effector catalogues of MRL8996, Fol4287 and Fo47 (CD-HIT at 70% identity) at the level of sequence, the candidates unique to MRL8996, their characterisation, and their structural comparison within the MRL8996 catalogue |
| 07 | Structural similarity network of the effector models including the unconnected nodes, so that the position of the unique candidates is visible |
| 08 | Maximum-likelihood phylogeny of the Fusarium oxysporum species complex |
| 09, 10, 19-22 | Per-chromosome content, repeat classes in the core and accessory compartments, effector genomic map and the composed effectorome figure |
| 11, 13-18 | Hi-C contact map, StainedGlass self-identity and GC criterion used to call the candidate centromeric regions, and the circos and Hi-C figure panels |
| 12, 23, 24 | Manuscript and response-letter documents |
| 25, 26 | Assembly of the data deposition (Data File 1-12 and all source data) and the rebuilt supplementary-table workbook |

The twelve large data tables previously distributed as manuscript files are
deposited as Data File 1-12 in the figshare record of this study (DOI to be
inserted on acceptance), together with the source data of every figure and
table.
