#!/usr/bin/env bash
# Stage 01.2 - Hi-C scaffolding (Juicer + 3D-DNA) to chromosome scale.
# Draft contigs + Hi-C reads -> scaffolds -> manual Juicebox review -> final AGP/FASTA.
# This documents the command flow; fill in the paths to your Juicer / 3D-DNA installs.
set -euo pipefail

DRAFT="FoxMRL8996_AssemblyScaffolds.fasta"
HIC_R1="[FILL hic_R1.fastq.gz]"
HIC_R2="[FILL hic_R2.fastq.gz]"
ENZYME="[FILL DpnII | Arima | none]"     # Hi-C restriction enzyme used in library prep
THREADS=32

# 1. reference prep
bwa index "$DRAFT"
samtools faidx "$DRAFT"
# python /path/to/juicer/misc/generate_site_positions.py "$ENZYME" MRL8996 "$DRAFT"

# 2. Juicer: map Hi-C reads -> aligned/merged_nodups.txt
#    (expects references/, restriction_sites/, fastq/ with HIC_R1/HIC_R2)
# juicer.sh -g MRL8996 -d "$PWD" -s "$ENZYME" -z references/"$DRAFT" \
#   -p references/"$DRAFT".chrom.sizes -y restriction_sites/MRL8996_"$ENZYME".txt \
#   -D /path/to/juicer -t "$THREADS"

# 3. 3D-DNA scaffolding -> *.assembly + *.hic  (your mrl8996-hic.assembly / mrl8996-hic.hic)
# run-asm-pipeline.sh "$DRAFT" aligned/merged_nodups.txt

# 4. after manual review in Juicebox Assembly Tools:
# run-asm-pipeline-post-review.sh -r mrl8996-hic.final.review.assembly \
#   "$DRAFT" aligned/merged_nodups.txt
#   -> mrl8996-hic.final.fasta / mrl8996-hic.final.agp

# 5. finalize chromosome naming/order + gap coordinates (your helper scripts)
# python format_final_fasta.py   # -> MRL8996_chromosome_level.fasta
# python find_gaps.py            # -> gaps_coordinates.bed
echo "Fill in Juicer/3D-DNA paths above, then run stepwise."
