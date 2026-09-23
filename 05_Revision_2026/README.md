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
| 01-04, 07, 09-12 | `python` (base conda) | pandas, matplotlib, numpy, biopython; CD-HIT and TM-align on `PATH` for 01 and 04 |
| 08 | `chr10-phylo` | blast, mafft, trimal, iqtree, pandas, matplotlib, biopython |
| 11a | `python` | `hicstraw` (pip), reads `mrl8996-hic.hic` directly |

Activate with `conda activate chr10-phylo` for script 08, `conda activate
base` for the rest.

## Run order and status

| # | Script | Produces | Status |
|---|---|---|---|
| 01 | `01_pan_effectorome_cdhit.py` | `results/pan_effectorome_MRL8996.fasta`, `clustered_effectors_MRL8996.clstr`, `tables/cdhit_cluster_matrix.tsv`, `tables/venn_counts.tsv`, `results/lists/*.txt`, `logs/cdhit.log` | done |
| 02 | `02_venn_effectorome.py` | `figures/Venn_Effectorome_MRL8996.{pdf,png}` | done |
| 03 | `03_characterise_unique_effectors.py` | `tables/Table_MRL8996_unique_effectors.tsv`, `tables/unique_effector_enrichment.tsv`, `figures/Fig_unique_effector_distribution.{pdf,png}` | done |
| 04 | `04_structural_homologs_within_MRL8996.py` | `tables/unique_vs_MRL_structural_homologs.tsv`, `unique_vs_MRL_structural_summary.tsv`, `figures/Fig_unique_intra_structural.{pdf,png}` | done |
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
- `08_FOSC_phylogeny.py --threads 8` sets BLAST/MAFFT/IQ-TREE threads.
- `11b_reexport_main_figures.py --dpi 600` sets the raster resolution;
  default is 600.
- `08_FOSC_phylogeny.py --plot-only` redraws the tree from the existing
  `FOSC_ML.treefile` without re-running BLAST/MAFFT/IQ-TREE.
- `12_build_manuscript_docx.py` is the one to re-run after you edit
  `manuscript/ms_v2.0_revised.md` or `response_to_comments.md`; it
  rebuilds both Word files and re-reads the tables from `tables/`.


## Scripts 13-26

### Centromeres (13-17)

| # | Script | Env | Output |
|---|--------|-----|--------|
| 13 | `13_stainedglass_percontig.py` | `sglass` | per-chromosome StainedGlass self-identity, `results/stainedglass/*.bed.gz` (16 files) |
| 14 | `14_call_centromeres_stainedglass.py` | `sglass` | first calls from the self-identity signal, compared with the Hi-C blocks of script 11a: `tables/Table_centromeres_stainedglass.tsv`, `results/centromeres_stainedglass.bed`, `figures/Fig_stainedglass_identity_16chr.{pdf,png}`. Superseded by script 17; kept for comparison |
| 17 | `17_centromeres_gc_criterion.py` | `python` | final calls: GC trough of each chromosome supported by a diverged identity band: `tables/Table_centromeres_gc.tsv`, `results/centromeres_gc.bed`, `tables/Table_identity_bands_per_chromosome.tsv`, `figures/Fig_centromere_identity_gc_profiles.{pdf,png}` |
| 15 | `15_circos_with_centromeres.R` | `r` (circlize) | circos plot with the script-17 centromeres: `figures/Figure_1B_circos_centromeres.pdf` |
| 16 | `16_label_hic_panelA.py` | `python` (pymupdf) | the original Juicebox export of the Hi-C map, cropped, with Chr01-Chr16 labels: `figures/Figure_1A_HiC_labelled.pdf` |

Run order: 11a -> 13 -> 14 -> 17 -> 15 -> 16 -> 11b.
Figure 1A is the labelled original Juicebox export from script 16; the map
re-rendered by script 11a is no longer used in Figure 1, but 11a still
provides `tables/Table_centromeric_interaction_blocks.tsv` (cis/trans Hi-C
signals per chromosome), which is Supplementary Table 1.

### Figures and tables (18-22)

| # | Script | Env | Output |
|---|--------|-----|--------|
| 18 | `18_supp_image1_tapestry_legend.py` | `python` | Supplementary Image 1 with read-depth key, telomere key, chromosome names and axis title; depth from `samtools bedcov` on the tapestry BAM: `figures/Supplementary_Image_1_tapestry.{pdf,png}`, `tables/Table_S5_read_depth_10kb.tsv` |
| 19 | `19_te_table_integrated_classes.py` | `python` | Table 2 on the integrated EDTA + DeepTE classification, core vs accessory: `tables/Table_TE_classes_core_vs_accessory.tsv` |
| 20 | `20_effector_map_plus_venn.py` | `python` | the authors' linear effector map (placed as vector) with the three-strain Venn diagram: `figures/Figure_effector_map_and_venn.{pdf,png}` |
| 21 | `21_effector_map_unique.R` | `r` | the authors' linear effector map with the 122 modelled MRL8996-unique effectors marked: `figures/MRL8996_Linear_Effector_Map_unique.pdf` |
| 22 | `22_figure4_recompose_network.py` | `python` | Figure 4 with panels A and C from `Figure_1_Master_Effectorome.pdf` unchanged and panel B replaced by the network of script 07 (`--panel`): `figures/Figure_4_Effectorome_recomposed.{pdf,png}` |

### Manuscript, letter and data deposition (23-26)

| # | Script | Env | Output |
|---|--------|-----|--------|
| 23 | `23_patch_v20_highlighted.py` | `python` | the authors' v2.0 manuscript with the revision edits in orange: `manuscript/MRL8996_2026_v2.1_highlighted.docx` |
| 24 | `24_build_response_letter.py` | `python` | point-by-point response letter: `manuscript/Response_to_reviewers_MRL8996_SciData.docx` |
| 25 | `25_build_deposition_package.py` | `python` | `deposition/`: Data File 1-12 and all source data, with `README.md` and a checksummed `MANIFEST.tsv` |
| 26 | `26_build_supplementary_workbook.py` | `python` | `tables/Supplementary_Tables_MRL8996_revised.xlsx` (Supplementary Tables 1-3) |

## Withdrawn analyses

The structural comparison of the isolate-specific candidates against the
Fol4287 and Fo47 catalogues (former scripts 05 and 06 and their outputs) has
been withdrawn: no structural models were produced for those two strains as
part of this study. The files are kept, unused, in
`withdrawn_Fo47_Fol4287_structural/`. The structural comparison reported in the
manuscript is the one within the MRL8996 catalogue (script 04) and the
three-strain comparison is at the level of sequence only (scripts 01-03).
