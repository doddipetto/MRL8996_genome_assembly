# *Fusarium oxysporum* MRL8996 — chromosome-level genome assembly

Code and annotation files accompanying the data descriptor
"A near-complete genome assembly of the *Fusarium oxysporum* keratitis isolate MRL8996".

The genome assembly and raw reads are deposited at NCBI (BioProject `PRJNA554890`):

- Genome assembly: GenBank `GCA_009746015.2`
- Nanopore reads: SRA `SRR39391754`
- Hi-C reads: SRA `SRR39391753`

A copy of the assembly FASTA is included here because every coordinate in the files below refers to it.

## Repository structure

```
01_Genome_Assembly/
  01_flye_assembly.sh              ONT de novo assembly (Flye)
  02_hic_scaffolding_juicer.sh     Hi-C scaffolding (Juicer + 3D-DNA)
  03_edta_te_annotation.sh         transposable-element annotation (EDTA)
02_Figure1_and_centromeres/
  15_circos_with_centromeres.R     Figure 1B (genome-feature circos plot)
  16_label_hic_panelA.py           Figure 1A (Hi-C contact map labels)
  17_centromeres_gc_criterion.py   candidate centromeric regions (StainedGlass + GC minima)
data/
  01_assembly_and_annotation/      assembly FASTA, AGP, gaps, lifted-over gene annotation (GFF3, BED),
                                   proteome, assembly statistics, per-chromosome content, 20-kb feature tracks
  02_repeats/                      EDTA TE annotation, TE library, DeepTE classification, repeat BED files,
                                   repeat content by core/accessory compartment
  03_centromeres_telomeres/        candidate centromeric regions and supporting signals; telomere calls (Tapestry)
  04_biosynthetic_gene_clusters/   antiSMASH region summary and region files
  05_technical_validation/         BUSCO, Merqury and QUAST summaries, read depth, auN, Hi-C contact matrix (.hic)
```

Tables are tab-separated; large files are gzip-compressed. Annotations are GFF3 or BED on the chromosome-level coordinates.

## Software

Flye, BWA, samtools, Juicer, 3D-DNA, Juicebox, TGS-GapCloser, liftoff, EDTA, RepeatMasker, DeepTE,
antiSMASH, StainedGlass, Tapestry, BUSCO, Merqury, QUAST; R (circlize). Versions and non-default
parameters are given in the Methods of the article; default parameters were used otherwise.

## Citation

Doddi, A., Puebla Planas, G., López-Berges, M. S. & Di Pietro, A. (2026). A near-complete genome assembly of the
*Fusarium oxysporum* keratitis isolate MRL8996.
