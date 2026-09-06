# Revision_2026 scripts

Everything for the revision lives under
`MRL8996_genome_paper/Revision_2026/` and writes only inside that folder
(`tables/`, `results/`, `figures/`, `data/`, `logs/`, `manuscript/`).
No script touches the original project directories — all inputs are read
read-only from `MRL8996_HiC_data/`, `MRL8996_genome_paper/` and
`F012_HiC_data/`.

Every script is standalone (`python3 NN_name.py`), takes no required
arguments, and re-runs idempotently: it overwrites its own outputs in
`Revision_2026/` and nothing else. Paths are hard-coded at the top of
each file under a `# --- paths ---` block if you need to move things.

## Environments

| Scripts | Environment | Needs |
|---|---|---|
| 01-07, 09-12 | `python` (base conda) | pandas, matplotlib, numpy, biopython; CD-HIT and TM-align on `PATH` for 01 and 04/05 |
| 08 | `chr10-phylo` | blast, mafft, trimal, iqtree, pandas, matplotlib, biopython |
| 11a | `python` | `hicstraw` (pip), reads `mrl8996-hic.hic` directly |
| `run_colabfold_Fo47.sh` | local GPU, ColabFold install | not yet run |

Activate with `conda activate chr10-phylo` for script 08, `conda activate
base` for the rest.

## Run order and status

| # | Script | Produces | Status |
|---|---|---|---|
| 01 | `01_pan_effectorome_cdhit.py` | `results/pan_effectorome_MRL8996.fasta`, `clustered_effectors_MRL8996.clstr`, `tables/cdhit_cluster_matrix.tsv`, `tables/venn_counts.tsv`, `results/lists/*.txt`, `logs/cdhit.log` | done |
| 02 | `02_venn_effectorome.py` | `figures/Venn_Effectorome_MRL8996.{pdf,png}` | done |
| 03 | `03_characterise_unique_effectors.py` | `tables/Table_MRL8996_unique_effectors.tsv`, `tables/unique_effector_enrichment.tsv`, `figures/Fig_unique_effector_distribution.{pdf,png}` | done |
| 04 | `04_structural_homologs_within_MRL8996.py` | `tables/unique_vs_MRL_structural_homologs.tsv`, `unique_vs_MRL_structural_summary.tsv`, `figures/Fig_unique_intra_structural.{pdf,png}` | done |
| 05 | `05_cross_strain_structural_Fol4287.py` | `results/tmalign_MRL8996unique_vs_Fol4287.tsv`, `tables/cross_strain_best_hits.tsv`, `figures/Fig_cross_strain_structural.{pdf,png}` | done |
| 06 | `06_prepare_Fo47_folding.py` + `run_colabfold_Fo47.sh` | `data/Fo47_effectors_for_colabfold.fasta`, `data/Fo47_header_map.tsv` | **fasta prepared, folding NOT run** |
| 07 | `07_structural_network_with_singletons.py` (`--panel` for the Figure 4 version) | `tables/network_singletons_nodes.tsv`, `network_family_composition.tsv`, `figures/Fig_structural_network_singletons.{pdf,png}` | done |
| 08 | `08_FOSC_phylogeny.py` | `results/phylogeny/FOSC_ML.*`, `tables/Table_FOSC_phylogeny_taxa.tsv`, `figures/Figure_FOSC_phylogeny.{pdf,png}`, `logs/08_phylogeny.log` | done |
| 09 | `09_accessory_chromosome_content.py` | `tables/Table_per_chromosome_content.tsv`, `Table_TE_classes_core_vs_accessory.tsv`, `effector_density_by_compartment.tsv`, `assembly_auN.tsv`, `figures/Fig_core_vs_accessory_content.{pdf,png}` | done |
| 10 | `10_effector_genomic_map.py` | `figures/Figure_effector_genomic_distribution.{pdf,png}` | done |
| 11a | `11a_hic_contact_map_centromeres.py` | `figures/Figure_1A_HiC_contact_map.{pdf,png}`, `tables/Table_centromeric_interaction_blocks.tsv` | done |
| 11b | `11b_reexport_main_figures.py` | `figures/final/Figure_{1..4}_*.{pdf,png}`, `figure_panel_provenance.tsv` | done |
| 12 | `12_build_manuscript_docx.py` | `manuscript/MRL8996_2026_v2.0_revised.docx`, `Response_to_comments.docx`, `tables/Table1_assembly_sequencing.tsv`, `Table5_three_strain.tsv` | done |

## Useful flags

- `11a_hic_contact_map_centromeres.py --binsize 50000` sets the contact-map
  bin size (default 50 kb). The `.hic` extraction is the slow part of this
  script; a coarser bin size runs faster.
- `05_cross_strain_structural_Fol4287.py --threads 8 --plot-only`: `--plot-only`
  re-summarises the existing TM-align pair table without repeating the
  61,000 alignments; `--threads` sets the alignment parallelism.
- `08_FOSC_phylogeny.py --threads 8` sets BLAST/MAFFT/IQ-TREE threads.
- `11b_reexport_main_figures.py --dpi 600` sets the raster resolution;
  default is 600.
- `08_FOSC_phylogeny.py --plot-only` redraws the tree from the existing
  `FOSC_ML.treefile` without re-running BLAST/MAFFT/IQ-TREE.
- `12_build_manuscript_docx.py` is the one to re-run after you edit
  `manuscript/ms_v2.0_revised.md` or `response_to_comments.md`; it
  rebuilds both Word files and re-reads the tables from `tables/`.

## The one thing still pending: Fo47 structural models

Script 06 wrote `data/Fo47_effectors_for_colabfold.fasta` (the Fo47
effector candidates with no existing model) and a header map. To fold
them on the local GPU:

```bash
cd /home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026
bash scripts/run_colabfold_Fo47.sh
``` Models are written to
`results/colabfold_Fo47/` (3 models per sequence, rank_001 retained), and
`COLABFOLD_BIN` at the top of the script currently points at
`/home/usuario/localcolabfold/.pixi/envs/default/bin/colabfold_batch` —
check that path before running.

Extending the cross-strain comparison to Fo47 afterwards needs a small
edit to script 05, not a flag: add the Fo47 model directory to its
`REF_SETS` block at the top and re-run it. The manuscript currently states in Methods that the
comparison is restricted to Fol4287 because no equivalent Fo47 models
were available, so that sentence would need updating too.

This run is not required for the revision as submitted — every figure,
table and number in `ms_v2.0_revised.md` comes from the scripts marked
"done" above.

## 14-15 (centromeres, re-run 2026-09-03)

| # | Script | Env | Status |
|---|--------|-----|--------|
| 13 | `13_stainedglass_percontig.py` | `sglass` | done - per-chromosome StainedGlass, 16 bed.gz in `results/stainedglass/` |
| 14 | `14_call_centromeres_stainedglass.py` | `sglass` | done - `tables/Table_centromeres_stainedglass.tsv`, `results/centromeres_stainedglass.bed`, `figures/Fig_stainedglass_identity_16chr.pdf/.png` |
| 15 | `15_circos_with_centromeres.R` | `r` (circlize) | done - `figures/Figure_1B_circos_centromeres.pdf` |

| 16 | `16_label_hic_panelA.py` | `python` (pymupdf) | done - `figures/Figure_1A_HiC_labelled.pdf`, Chr01-Chr16 measured from the blue boxes |
| 17 | `17_centromeres_gc_criterion.py` | `python` | done - final centromere calls from GC minima with identity support; supersedes script 14. Writes `tables/Table_centromeres_gc.tsv`, `results/centromeres_gc.bed`, `tables/Table_identity_bands_per_chromosome.tsv`, `figures/Fig_centromere_identity_gc_profiles.*`. Run before 15 and 11b |
| 18 | `18_supp_image1_tapestry_legend.py` | `python` | done - adds the read-depth key, telomere key, chromosome names and axis title to the tapestry contig plot (R2.21); depth analysis behind the legend text is in results/coverage (samtools bedcov on the tapestry BAM) -> figures/Supplementary_Image_1_tapestry.{png,pdf}, tables/Table_S5_* |
| 19 | `19_te_table_integrated_classes.py` | `python` | done - rebuilds Table 2 on the integrated EDTA + DeepTE classification, pooled to class/superfamily, core vs accessory bp and %, subtotals and non-redundant totals in-table (R2.20) -> tables/Table_TE_classes_core_vs_accessory.tsv |
| 20 | `20_effector_map_plus_venn.py` | `python` | done - two-panel figure: author's linear effector map (panel a, placed as vector, not redrawn) + CD-HIT pan-effectorome Venn (panel b) -> figures/Figure_effector_map_and_venn.pdf/.png |
| 21 | `21_effector_map_unique.R` | `r` | done - author's linear effector map with the 122 modelled MRL8996-unique effectors marked (black triangles); geometry/theme copied verbatim from plot_map_effector_position.R -> figures/MRL8996_Linear_Effector_Map_unique.pdf |
| 22 | `22_figure4_recompose_network.py` | `python` | done - Figure 4 rebuilt: panels A and C clipped from `Figure_1_Master_Effectorome.pdf` unchanged, panel B replaced by the singleton-retaining network (`07 --panel`) -> `figures/Figure_4_Effectorome_recomposed.{pdf,png}`; `11b` now re-exports Figure 4 from it |
| 23 | `23_patch_v20_highlighted.py` | patches the author v2.0 docx in place with the corrections still missing, marked in orange | `MRL8996_2026_v2.1_highlighted.docx` |
| 24 | `24_build_response_letter.py` | point-by-point response letter in the FO12 Scientific Data style | `Response_to_reviewers_MRL8996_SciData.docx` |

Run order: 13 -> 14 -> 15 -> 16 -> `11b_reexport_main_figures.py` -> `12_build_manuscript_docx.py`.
Figure 1 panel A is now the cropped original Juicebox export with chromosome
labels added by script 16 (`figures/Figure_1A_HiC_labelled.pdf`) rather than the re-rendered map
from script 11a; script 11a is retained but no longer feeds Figure 1.
Table 3 in the docx builder now reads `Table_centromeres_stainedglass.tsv`; the
old Hi-C table becomes Supplementary Table S4.
| 25 | `25_build_deposition_package.py` | assembles `deposition/` (the twelve large data tables as Data File 1-12 plus all source data), writes `README.md` and a checksummed `MANIFEST.tsv` | `deposition/` |
| 26 | `26_build_supplementary_workbook.py` | rebuilds the supplementary-table workbook with only the three small tables that stay with the manuscript | `tables/Supplementary_Tables_MRL8996_revised.xlsx` |
