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

- **Mature sequences for structure prediction.** SignalP-cleaved (mature) sequences
  are folded in `02_Effector_Prediction/04`. AMAPEC explicitly requires structures
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

Flye, Medaka (optional), BWA, samtools, Juicer, 3D-DNA, Juicebox; EDTA (bedtools); SignalP 6,
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
