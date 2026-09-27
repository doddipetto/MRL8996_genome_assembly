# 05_Revision_2026

Material added during the revision of the MRL8996 genome paper.

## Data_tables/

The data tables as Excel files: Data File 1-12 (the tables cited in the
manuscript as Data File 1 to 12) and `Structural_network_MRL8996.xlsx`
(the structural network of Figure 5B, with the closest known fold of each
family and the proteins unique to MRL8996). See `Data_tables/README.md`.

## Scripts for the updated figures

| Figure | Script | Output |
|---|---|---|
| Figure 1A | `16_label_hic_panelA.py` | Hi-C contact map (original Juicebox export) with Chr01-Chr16 labels |
| Figure 1B | `15_circos_with_centromeres.R` | circos plot with the candidate centromeric regions |
| Figure 4A | `21_effector_map_unique.R` | linear effector map with the MRL8996-unique candidates marked |
| Figure 4B | `02_venn_effectorome.py` | three-way Venn diagram of the MRL8996, Fol4287 and Fo47 effector catalogues (CD-HIT, 70% identity) |
| Figure 4 | `20_effector_map_plus_venn.py` | assembles panels A and B |
| Figure 5B | `07_structural_network_with_singletons.py --panel` | structural similarity network with all 482 nodes, including the singletons |
| Figure 5 | `22_figure4_recompose_network.py` | replaces panel B of the effectorome figure with the new network |
| Supplementary Image 1 | `18_supp_image1_tapestry_legend.py` | tapestry telomere plot with read-depth and telomere keys |
| Supplementary Image 4 | `17_centromeres_gc_criterion.py` | per-chromosome GC and self-identity profiles; also writes the centromere calls read by `15` |

Run 17 before 15, 21 and 02 before 20, and 07 before 22. Paths inside the
scripts are absolute to the analysis machine and have to be adapted. Tools:
python (pandas, matplotlib, numpy, networkx, pymupdf, reportlab), R (circlize,
ggplot2).
