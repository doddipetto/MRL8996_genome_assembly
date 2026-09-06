#!/usr/bin/env python3
"""Assemble the repository deposition for the MRL8996 data descriptor.

Editor comment E7 asks for the large data tables to be taken out of the
manuscript files and deposited in a repository. This script collects them,
together with the source data a reader needs to reproduce every figure and
table, into a single directory tree with a README and a checksummed manifest.

Nothing outside the deposition directory is written; all sources are copied,
never moved.

out: Revision_2026/deposition/            tree ready to upload
     Revision_2026/deposition/README.md   dataset description
     Revision_2026/deposition/MANIFEST.tsv  file, bytes, md5, description
"""
import csv
import gzip
import hashlib
import shutil
import subprocess
import tarfile
from pathlib import Path

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
HIC = ADODDI / "MRL8996_HiC_data"
PAPER = ADODDI / "MRL8996_genome_paper"
REV = PAPER / "Revision_2026"
JGI = (PAPER / "JGI_reference" / "FoxMRL8996" / "Mycocosm" / "Annotation"
       / "Filtered_Models___best__")
PROT = JGI / "Proteins"
OUT = REV / "deposition"
SUPP_XLSX = PAPER / "media-1 (3).xlsx"   # the submitted supplementary-table workbook

# (destination relative to OUT, source, gzip?, description)
# Destinations named Data_File_NN_* are the twelve tables previously supplied as
# manuscript files; the numbering follows the in-text citations of the revised
# manuscript.
FILES = [
    # ------------------------------------------------ 01 assembly + annotation
    ("01_assembly_and_annotation/assembly_and_sequencing_statistics.tsv",
     REV / "tables" / "Table1_assembly_sequencing.tsv", False,
     "Sequencing data and assembly statistics of MRL8996 v1.0 and "
     "v2.0, including contig and scaffold N50/L50, auN, ambiguous bases, BUSCO "
     "and Merqury values."),
    ("01_assembly_and_annotation/MRL8996_chromosome_level.fasta.gz",
     HIC / "MRL8996_chromosome_level.fasta", True,
     "Chromosome-level assembly (16 pseudomolecules) as used for every analysis "
     "in the manuscript; identical to GenBank GCA_009746015.2."),
    ("01_assembly_and_annotation/MRL8996_chromosome_level.fasta.fai",
     HIC / "MRL8996_chromosome_level.fasta.fai", False,
     "Index of the assembly: pseudomolecule names and lengths."),
    ("01_assembly_and_annotation/mrl8996-hic.final.agp",
     HIC / "mrl8996-hic.final.agp", False,
     "AGP file describing how the Nanopore contigs were placed and oriented "
     "into the 16 pseudomolecules by Hi-C scaffolding."),
    ("01_assembly_and_annotation/gaps_coordinates.bed",
     HIC / "gaps_coordinates.bed", False,
     "Coordinates of the gap-spanning N runs introduced during scaffolding."),
    ("01_assembly_and_annotation/MRL8996_chromosome_level.gff3.gz",
     HIC / "08_circos_features" / "MRL8996_chromosome_level.gff3", True,
     "Gene annotation in GFF3 format on the chromosome-level coordinates, "
     "lifted over from the v1.0 JGI filtered 'best' gene models with liftoff."),
    ("01_assembly_and_annotation/MRL8996_proteome.aa.fasta.gz",
     PROT / "FoxMRL8996_GeneCatalog_proteins_20200124.aa.fasta", True,
     "Predicted proteome of the filtered 'best' gene models (JGI Mycocosm "
     "annotation), the input of the effector pipeline."),
    ("01_assembly_and_annotation/MRL8996_genes.bed",
     HIC / "08_circos_features" / "MRL8996_genes.bed", False,
     "Gene positions in BED format, used for the gene-density track."),
    ("01_assembly_and_annotation/per_chromosome_content.tsv",
     REV / "tables" / "Table_per_chromosome_content.tsv", False,
     "Per-chromosome length, gene count and density, coding fraction, "
     "repeat-masked fraction, effector count and biosynthetic gene clusters, "
     "with the core or accessory assignment of each pseudomolecule."),
    ("01_assembly_and_annotation/mapping_cromosomi_NCBI.tsv",
     HIC / "mapping_cromosomi_NCBI.tsv", False,
     "Correspondence between the internal pseudomolecule names and the GenBank "
     "sequence accessions."),

    # ------------------------------------------------------------- 02 synteny
    ("02_synteny/synteny_matrix_MRL8996_vs_Fol4287.csv",
     HIC / "06_circos_synteny" / "Matrix_MRL_vs_Fol4287.csv", False,
     "Matrix of aligned sequence (bp) between each MRL8996 "
     "pseudomolecule and each Fol4287 chromosome."),
    ("02_synteny/synteny_matrix_MRL8996_vs_Fo47.csv",
     HIC / "06_circos_synteny" / "Matrix_Fo47_vs_MRL8996.csv", False,
     "Matrix of aligned sequence (bp) between each MRL8996 "
     "pseudomolecule and each Fo47 chromosome."),
    ("02_synteny/self_alignment_matrix.csv",
     HIC / "06_circos_synteny" / "Matrix_MRL_Self.csv", False,
     "Matrix of aligned sequence (bp) between MRL8996 "
     "pseudomolecules in the self-alignment (segmental duplications, "
     "Figure 3B)."),
    ("02_synteny/FO47_vs_MRL8996.paf",
     HIC / "06_circos_synteny" / "FO47_vs_MRL8996.paf", False,
     "minimap2 asm20 alignment of Fo47 against MRL8996, the source of Data "
     "File 3 and of the Figure 3A ribbons."),
    ("02_synteny/FO47_vs_FOL4287.paf",
     HIC / "06_circos_synteny" / "FO47_vs_FOL4287.paf", False,
     "minimap2 asm20 alignment of Fo47 against Fol4287."),
    ("02_synteny/synteny_JGI_vs_HiC.paf",
     HIC / "12_DGenies_Synteny" / "synteny_JGI_vs_HiC.paf", False,
     "minimap2 alignment of the v1.0 draft against the chromosome-level "
     "assembly, the source of the Figure 2A dot plot."),

    # ------------------------------------------------------------- 03 repeats
    ("03_repeats/TE_annotation_EDTA.gff3.gz",
     HIC / "09_TE_annotation" / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.TEanno.gff3",
     True,
     "Full transposable-element annotation produced by EDTA "
     "v2.2.2 on the chromosome-level assembly."),
    ("03_repeats/MRL8996_TE_library.fa.gz",
     HIC / "09_TE_annotation" / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.TElib.fa",
     True,
     "Non-redundant TE library built by EDTA and used to soft-mask the "
     "assembly with RepeatMasker."),
    ("03_repeats/MRL8996_TE_library_novel.fa.gz",
     HIC / "09_TE_annotation" / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.TElib.novel.fa",
     True,
     "Subset of the TE library not matching known Repbase families."),
    ("03_repeats/DeepTE_classification.txt.gz",
     HIC / "09_TE_annotation" / "DeepTE" / "DeepTE-master" / "output_dir"
     / "opt_DeepTE.txt", True,
     "DeepTE classification of the elements left unclassified by EDTA."),
    ("03_repeats/TE_classes_core_vs_accessory.tsv",
     REV / "tables" / "Table_TE_classes_core_vs_accessory.tsv", False,
     "Counts and occupancy (bp and per cent) of every TE class and superfamily "
     "in the core compartment, in the accessory compartment and genome-wide; "
     "the source of Table 2 of the manuscript."),
    ("03_repeats/MRL8996_complex_TEs.bed",
     HIC / "08_circos_features" / "MRL8996_complex_TEs.bed", False,
     "Positions of the complex (class I and class II) elements."),
    ("03_repeats/MRL8996_simple_repeats.bed",
     HIC / "08_circos_features" / "MRL8996_simple_repeats.bed", False,
     "Positions of the simple and low-complexity repeats."),

    # ----------------------------------------------------------- 04 effectors
    ("04_effectors/structural_families_per_candidate.tsv",
     PROT / "Final_Effectorome_Families_Complex.tsv", False,
     "Per-candidate table of the structural classification: "
     "family assignment from the TM-score dendrogram, cluster size, pLDDT, "
     "AMAPEC prediction and annotation of the 482 modelled candidates."),
    ("04_effectors/structural_network_nodes.tsv",
     PAPER / "structural_network" / "OCE_struct_network_nodes.tsv", False,
     "Node table of the force-directed structural-similarity "
     "network: candidate, community, degree and layout coordinates."),
    ("04_effectors/structural_network_nodes_with_singletons.tsv",
     REV / "tables" / "network_singletons_nodes.tsv", False,
     "Node table of the revised network in which all 482 nodes are retained, "
     "including the 131 with no edge above TM-score 0.5 (Figure 5B)."),
    ("04_effectors/structural_network_family_composition.tsv",
     REV / "tables" / "network_family_composition.tsv", False,
     "Composition of each structural family in the revised network, with the "
     "number of isolate-specific candidates per family."),
    ("04_effectors/structural_network.gml",
     PAPER / "structural_network" / "OCE_struct_network.gml", False,
     "The structural-similarity network in GML format, readable in Cytoscape "
     "or Gephi."),
    ("04_effectors/SignalP_whole_proteome.csv.gz",
     PROT / "SignalP_results" / "MRL8996_SignalP_Table.csv", True,
     "SignalP v6.0 predictions for the whole predicted proteome, "
     "with the predicted signal-peptide class and cleavage site."),
    ("04_effectors/EffectorP_whole_proteome.csv.gz",
     PROT / "MRL8996_EffectorP_Table.csv", True,
     "EffectorP v3.0 predictions for the whole predicted "
     "proteome, with the apoplastic, cytoplasmic and non-effector "
     "probabilities."),
    ("04_effectors/SignalP_signal_peptide_regions.gff3.gz",
     PROT / "SignalP_results" / "region_output.gff3", True,
     "Signal-peptide regions predicted by SignalP in GFF3 "
     "format."),
    ("04_effectors/TM_score_all_vs_all_full.tsv.gz",
     PROT / "tm_align_results_all_vs_all.tsv", True,
     "All-versus-all TM-align comparison of the 482 modelled "
     "candidates (Protein_A, Protein_B, Max_TM_Score), the input of both the "
     "heatmap and the network."),
    ("04_effectors/effector_catalogue_annotated.tsv",
     HIC / "effector_map_data.tsv", False,
     "The 513 candidate effectors with their genomic position, "
     "compartment, model confidence (pLDDT), AMAPEC prediction and "
     "antimicrobial probability; the source of Figure 4A and of the AMAPEC "
     "track of Figure 5A."),
    ("04_effectors/MRL8996_final_effectorome.fasta",
     PROT / "FoxMRL8996_Final_Effectorome.fasta", False,
     "Amino-acid sequences of the 513 candidate effectors."),
    ("04_effectors/MRL8996_SignalP_EffectorP_intersect.fasta",
     PROT / "SignalP_EffectorP_intersect.fasta", False,
     "Sequences shared by the SignalP and EffectorP predictions, the starting "
     "set before the DeepTMHMM filter."),
    ("04_effectors/MRL8996_final_effectorome_IDs.txt",
     PROT / "FoxMRL8996_Final_Effectorome_IDs.txt", False,
     "JGI protein identifiers of the 513 candidates."),
    ("04_effectors/effector_density_by_compartment.tsv",
     REV / "tables" / "effector_density_by_compartment.tsv", False,
     "Candidate counts and densities in the core and accessory compartments."),

    # ------------------------------------------- 05 three-strain comparison
    ("05_three_strain_comparison/cdhit_clusters.clstr",
     REV / "results" / "clustered_effectors_MRL8996.clstr", False,
     "CD-HIT clustering (70% identity) of the pooled effector catalogues of "
     "MRL8996, Fol4287 and Fo47, in native .clstr format."),
    ("05_three_strain_comparison/cdhit_cluster_matrix.tsv",
     REV / "tables" / "cdhit_cluster_matrix.tsv", False,
     "The same clustering as a table: one row per sequence with its cluster, "
     "strain and representative status."),
    ("05_three_strain_comparison/pan_effectorome.fasta",
     REV / "results" / "pan_effectorome_MRL8996.fasta", False,
     "The 1,406 pooled effector sequences of the three strains."),
    ("05_three_strain_comparison/venn_counts.tsv",
     REV / "tables" / "venn_counts.tsv", False,
     "Cluster and protein counts of the seven regions of the three-way "
     "comparison (Figure 4B)."),
    ("05_three_strain_comparison/MRL8996_unique_effectors.tsv",
     REV / "tables" / "Table_MRL8996_unique_effectors.tsv", False,
     "The 135 isolate-specific candidates with cluster, chromosome, "
     "compartment, pLDDT, AMAPEC prediction and structural family."),
    ("05_three_strain_comparison/unique_effector_enrichment.tsv",
     REV / "tables" / "unique_effector_enrichment.tsv", False,
     "Compartment and antimicrobial-prediction enrichment tests for the "
     "isolate-specific set."),
    ("05_three_strain_comparison/unique_vs_MRL8996_structural_homologs.tsv",
     REV / "tables" / "unique_vs_MRL_structural_homologs.tsv", False,
     "Best structural counterpart of each isolate-specific candidate within "
     "the MRL8996 catalogue, with the TM-score."),
    ("05_three_strain_comparison/unique_vs_MRL8996_structural_summary.tsv",
     REV / "tables" / "unique_vs_MRL_structural_summary.tsv", False,
     "Summary of the same comparison at several TM-score thresholds."),
    ("05_three_strain_comparison/unique_vs_Fol4287_TM_scores.tsv.gz",
     REV / "results" / "tmalign_MRL8996unique_vs_Fol4287.tsv", True,
     "All TM-align comparisons between the isolate-specific candidates and the "
     "Fol4287 models."),
    ("05_three_strain_comparison/cross_strain_best_hits.tsv",
     REV / "tables" / "cross_strain_best_hits.tsv", False,
     "Best cross-strain structural hit of every MRL8996 candidate."),
    ("05_three_strain_comparison/three_strain_summary.tsv",
     REV / "tables" / "Table5_three_strain.tsv", False,
     "Catalogue sizes and clustering summary of the three strains."),
    ("05_three_strain_comparison/Fo47_effectors_for_colabfold.fasta",
     REV / "data" / "Fo47_effectors_for_colabfold.fasta", False,
     "Fo47 effector sequences submitted to ColabFold for the cross-strain "
     "structural comparison."),
    ("05_three_strain_comparison/Fo47_header_map.tsv",
     REV / "data" / "Fo47_header_map.tsv", False,
     "Mapping between the shortened ColabFold headers and the original Fo47 "
     "accessions."),

    # --------------------------------------------------------- 06 centromeres
    ("06_centromeres/centromeric_interaction_blocks.tsv",
     REV / "tables" / "Table_centromeric_interaction_blocks.tsv", False,
     "Candidate centromeric region of each pseudomolecule with the supporting "
     "signals (GC minimum, Hi-C interaction focus, StainedGlass identity "
     "band)."),
    ("06_centromeres/centromeres_gc.tsv",
     REV / "tables" / "Table_centromeres_gc.tsv", False,
     "GC content in 20-kb windows and the position of the per-chromosome "
     "minimum."),
    ("06_centromeres/centromeres_stainedglass.tsv",
     REV / "tables" / "Table_centromeres_stainedglass.tsv", False,
     "Self-identity bands called from the StainedGlass maps."),
    ("06_centromeres/identity_bands_per_chromosome.tsv",
     REV / "tables" / "Table_identity_bands_per_chromosome.tsv", False,
     "Extent and identity range of the self-identity bands per chromosome."),
    ("06_centromeres/centromeres_gc.bed",
     REV / "results" / "centromeres_gc.bed", False,
     "GC-minimum intervals in BED format."),
    ("06_centromeres/centromeres_stainedglass.bed",
     REV / "results" / "centromeres_stainedglass.bed", False,
     "StainedGlass identity-band intervals in BED format."),

    # -------------------------------------------------- 07 technical validation
    ("07_technical_validation/BUSCO_short_summary.txt",
     HIC / "04_busco_analysis" / "MRL8996_BUSCO"
     / "short_summary.specific.hypocreales_odb10.MRL8996_BUSCO.txt", False,
     "BUSCO v5 summary against hypocreales_odb10."),
    ("07_technical_validation/merqury_qv.txt",
     HIC / "05_merqury_analysis" / "MRL8996_stats.qv", False,
     "Merqury consensus quality value."),
    ("07_technical_validation/merqury_completeness.stats",
     HIC / "05_merqury_analysis" / "MRL8996_stats.completeness.stats", False,
     "Merqury k-mer completeness."),
    ("07_technical_validation/quast_report.tsv",
     HIC / "assembly_comparison_quast" / "report.tsv", False,
     "QUAST comparison of the v1.0 draft and the chromosome-level assembly."),
    ("07_technical_validation/read_depth_10kb.tsv.gz",
     REV / "tables" / "Table_S5_read_depth_10kb.tsv", True,
     "Nanopore read depth in 10-kb windows, all alignments and MAPQ>=20, used "
     "to examine the regions of elevated coverage in Supplementary Image 1."),
    ("07_technical_validation/windows_above_2x.tsv",
     REV / "tables" / "Table_S5_windows_above_2x.tsv", False,
     "Windows exceeding twice the median depth, with their annotation (rDNA "
     "array, mitochondrial insertion, Chr13 first contig)."),
    ("07_technical_validation/assembly_auN.tsv",
     REV / "tables" / "assembly_auN.tsv", False,
     "auN of both assembly versions and the underlying length distribution."),
    ("07_technical_validation/mrl8996-hic.hic",
     HIC / "mrl8996-hic.hic", False,
     "Hi-C contact matrix in .hic format, readable in Juicebox; the source of "
     "Figure 1A."),

    # ------------------------------------------------------------ 08 phylogeny
    ("08_phylogeny/FOSC_concatenated_alignment.fasta.gz",
     REV / "results" / "phylogeny" / "concatenated.fasta", True,
     "Concatenated amino-acid alignment of the single-copy orthologues used "
     "for the FOSC maximum-likelihood phylogeny."),
    ("08_phylogeny/FOSC_ML.treefile",
     REV / "results" / "phylogeny" / "FOSC_ML.treefile", False,
     "Maximum-likelihood tree (IQ-TREE) of the FOSC strains including MRL8996, "
     "Fol4287 and Fo47."),
    ("08_phylogeny/FOSC_ML.contree",
     REV / "results" / "phylogeny" / "FOSC_ML.contree", False,
     "Consensus tree with ultrafast-bootstrap support."),
    ("08_phylogeny/FOSC_phylogeny_taxa.tsv",
     REV / "tables" / "Table_FOSC_phylogeny_taxa.tsv", False,
     "Strains included in the phylogeny with their assembly accessions."),
]


# The twelve tables submitted as a single workbook of manuscript files. They are
# exported sheet by sheet, verbatim, and become Data File 1-12 of the
# deposition; the manuscript now cites them under those names.
SHEETS = [
    (1, "additional_assembly_statistics", False,
     "Additional assembly statistics of MRL8996 (scaffold counts, length "
     "distribution, GC, N content)."),
    (2, "identity_matrix_MRL8996_vs_Fol4287", False,
     "Matrix of identity between MRL8996 and Fol4287, chromosome by "
     "chromosome."),
    (3, "identity_matrix_MRL8996_vs_Fo47", False,
     "Matrix of identity between MRL8996 and Fo47, chromosome by chromosome."),
    (4, "identity_matrix_MRL8996_self", False,
     "Matrix of identity of MRL8996 against itself (segmental duplications)."),
    (5, "structural_family_distribution", False,
     "Distribution of the structural protein families identified from the "
     "TM-score dendrogram, with the number of members per family."),
    (6, "structural_network_nodes", False,
     "Structural network nodes of the MRL8996 effector structures, with family, "
     "degree and layout coordinates."),
    (7, "repeat_statistics_per_chromosome", False,
     "Statistics of repeat sequences in the MRL8996 genome per chromosome."),
    (8, "SignalP_results", False,
     "SignalP v6.0 results for the predicted secreted proteins of MRL8996."),
    (9, "EffectorP_results", False,
     "EffectorP v3.0 results for the potential effector proteins of MRL8996."),
    (10, "DeepTMHMM_results", False,
     "DeepTMHMM results on the MRL8996 effector candidates (topology and "
     "presence of transmembrane domains)."),
    (11, "TM_score_all_vs_all", True,
     "All-versus-all TM-scores of the MRL8996 AlphaFold2 models "
     "(Protein_A, Protein_B, Max_TM_Score)."),
    (12, "AMAPEC_predictions", False,
     "Antimicrobial-activity prediction of the MRL8996 effector candidates "
     "with AMAPEC, with the per-model pLDDT."),
]

# Circos figure tracks, copied as a block.
TRACKS = [
    ("track_gc_content_20kb.txt", "GC content in 20-kb windows"),
    ("track_gene_density_20kb.txt", "Gene density in 20-kb windows"),
    ("track_exon_density_20kb.txt", "Exon density in 20-kb windows"),
    ("windows_20kb.bed", "The 20-kb window set used for all density tracks"),
    ("chrom.sizes", "Pseudomolecule lengths"),
]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def place(dest, src, do_gzip, rows, desc):
    src = Path(src)
    if not src.exists():
        print("  MISSING", src)
        return
    target = OUT / dest
    target.parent.mkdir(parents=True, exist_ok=True)
    if do_gzip:
        with open(src, "rb") as fi, gzip.open(target, "wb", compresslevel=6) as fo:
            shutil.copyfileobj(fi, fo, length=1 << 22)
    else:
        shutil.copy2(src, target)
    rows.append((dest, target.stat().st_size, md5(target), desc))


def export_supplementary(rows):
    """Export every sheet of the submitted workbook as a Data File."""
    import openpyxl
    if not SUPP_XLSX.exists():
        print("  MISSING", SUPP_XLSX)
        return
    wb = openpyxl.load_workbook(str(SUPP_XLSX), read_only=True, data_only=True)
    dest = OUT / "00_data_files_cited_in_manuscript"
    dest.mkdir(parents=True, exist_ok=True)
    for n, slug, do_gzip, desc in SHEETS:
        ws = wb["Supplementary Table %d" % n]
        name = "Data_File_%02d_%s.tsv" % (n, slug)
        target = dest / (name + ".gz" if do_gzip else name)
        op = gzip.open(target, "wt", newline="", compresslevel=6) if do_gzip \
            else open(target, "w", newline="")
        with op as fo:
            w = csv.writer(fo, delimiter="\t", lineterminator="\n")
            for r in ws.iter_rows(values_only=True):
                if all(c is None for c in r):
                    continue
                w.writerow(["" if c is None else c for c in r])
        rows.append((str(target.relative_to(OUT)), target.stat().st_size,
                     md5(target), "Data File %d. %s" % (n, desc)))
    wb.close()


def bgc_archive(rows):
    """Tar the antiSMASH region files and derive a summary table."""
    bgc = HIC / "10_BGC_annotation" / "MRL8996_BGC"
    gbks = sorted(bgc.glob("*.region*.gbk"))
    if not gbks:
        print("  MISSING", bgc)
        return
    dest = OUT / "09_biosynthetic_gene_clusters"
    dest.mkdir(parents=True, exist_ok=True)

    summary = dest / "BGC_regions_summary.tsv"
    with open(summary, "w") as fo:
        fo.write("chromosome\tregion\tstart\tend\tproduct\n")
        for g in gbks:
            chrom, region = g.name.split(".")[0], g.name.split(".")[1]
            start = end = ""
            products = []
            with open(g) as fh:
                for line in fh:
                    s = line.strip()
                    if s.startswith("/product=") and len(products) < 6:
                        products.append(s.split("=", 1)[1].strip('"'))
                    elif s.startswith("Orig. start"):
                        start = s.split("::")[1].strip()
                    elif s.startswith("Orig. end"):
                        end = s.split("::")[1].strip()
                    elif s.startswith("FEATURES") and start and end and products:
                        break
            fo.write("%s\t%s\t%s\t%s\t%s\n"
                     % (chrom, region, start, end, ";".join(dict.fromkeys(products))))
    rows.append((str(summary.relative_to(OUT)), summary.stat().st_size, md5(summary),
                 "Coordinates and predicted product of every biosynthetic gene "
                 "cluster detected by antiSMASH, one row per region."))

    tar = dest / "antiSMASH_regions_gbk.tar.gz"
    with tarfile.open(tar, "w:gz") as tf:
        for g in gbks:
            tf.add(g, arcname="antiSMASH_regions/" + g.name)
    rows.append((str(tar.relative_to(OUT)), tar.stat().st_size, md5(tar),
                 "Full antiSMASH output for every region in GenBank format "
                 "(%d files)." % len(gbks)))


README = """# Fusarium oxysporum MRL8996: chromosome-level assembly, annotation and \
effector catalogue - source data

This deposition accompanies the data descriptor "{title}". It holds the large
data tables cited in the manuscript, which are distributed here rather than as
manuscript files, together with the annotation and the intermediate results
needed to reproduce every figure and table.

The assembly itself and the raw reads are in the primary archives and are not
duplicated here:

- Genome assembly: GenBank GCA_009746015.2 (BioProject PRJNA554890)
- Nanopore reads: SRA SRR39391754
- Hi-C reads: SRA SRR39391753

A copy of the assembly FASTA is nevertheless included, because every coordinate
in the files below refers to it.

## Data files cited in the manuscript

The twelve tables previously supplied as manuscript files are named
`Data_File_NN_*` and are cited in the text as Data File 1 to Data File 12:

{datafiles}

## Directory contents

| Directory | Contents |
|---|---|
| `00_data_files_cited_in_manuscript` | Data File 1 to Data File 12, the twelve tables that were previously distributed as manuscript files, exported verbatim from the submitted workbook |
| `01_assembly_and_annotation` | Assembly FASTA and index, AGP and gap coordinates, gene annotation (GFF3) on the chromosome-level coordinates, predicted proteome, per-chromosome content table, accession mapping |
| `02_synteny` | Chromosome-by-chromosome alignment matrices against Fol4287 and Fo47 and for the self-alignment, and the underlying minimap2 PAF files |
| `03_repeats` | EDTA transposable-element annotation, non-redundant TE library, DeepTE classification, per-class occupancy in the core and accessory compartments, repeat BED files |
| `04_effectors` | Effector catalogue with all per-candidate scores, SignalP and EffectorP predictions for the whole proteome, all-versus-all TM-score matrix, structural family assignment, structural-similarity network (node table and GML) |
| `05_three_strain_comparison` | CD-HIT clustering of the MRL8996, Fol4287 and Fo47 catalogues, the seven-region counts of the three-way comparison, the isolate-specific candidates and their structural comparisons |
| `06_centromeres` | Candidate centromeric regions with their supporting signals, GC profiles, StainedGlass identity bands |
| `07_technical_validation` | BUSCO, Merqury and QUAST summaries, read depth in 10-kb windows, auN, Hi-C contact matrix (.hic) |
| `08_phylogeny` | Concatenated alignment and maximum-likelihood trees of the FOSC strains |
| `09_biosynthetic_gene_clusters` | antiSMASH region summary and the full region files |

`MANIFEST.tsv` lists every file with its size, MD5 checksum and a one-line
description.

## Formats

Tables are tab-separated unless the extension says otherwise; the synteny
matrices are comma-separated as produced. Files larger than a few megabytes are
gzip-compressed. Sequence files are FASTA, annotations are GFF3 or BED on the
chromosome-level coordinates, alignments are PAF, the network is GML and the
Hi-C matrix is in .hic format (Juicebox).

## Software that produced these files

The command lines are in the Methods of the manuscript and in the code
deposition. In short: Flye and Hi-C scaffolding for the assembly, liftoff for
the annotation transfer, EDTA and DeepTE for the repeats, SignalP, EffectorP,
DeepTMHMM, ColabFold (AlphaFold2), TM-align and AMAPEC for the effector
catalogue, CD-HIT for the cross-strain clustering, minimap2 for the alignments,
StainedGlass for the self-identity maps, IQ-TREE for the phylogeny and
antiSMASH for the biosynthetic gene clusters.
"""

TITLE = ("A near-complete genome assembly of the Fusarium oxysporum keratitis "
         "isolate MRL8996")


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    rows = []
    export_supplementary(rows)
    for dest, src, gz, desc in FILES:
        place(dest, src, gz, rows, desc)
    for name, desc in TRACKS:
        place("01_assembly_and_annotation/tracks/" + name,
              HIC / "08_circos_features" / name, False, rows, desc)
    bgc_archive(rows)

    datafiles = "\n".join(
        "- `%s` - %s" % (d, desc.split(". ", 1)[1])
        for d, _, _, desc in sorted(rows) if "Data_File_" in d)

    (OUT / "README.md").write_text(
        README.format(title=TITLE, datafiles=datafiles))

    with open(OUT / "MANIFEST.tsv", "w") as fo:
        fo.write("file\tbytes\tmd5\tdescription\n")
        for d, n, h, desc in sorted(rows):
            fo.write("%s\t%d\t%s\t%s\n" % (d, n, h, desc))

    total = sum(n for _, n, _, _ in rows)
    print("files: %d | total: %.1f MB | %s"
          % (len(rows), total / 1e6, OUT))
    print(subprocess.run(["du", "-sh", str(OUT)], capture_output=True,
                         text=True).stdout.strip())


if __name__ == "__main__":
    main()
