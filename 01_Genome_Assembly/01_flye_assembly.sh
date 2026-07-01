#!/usr/bin/env bash
# Stage 01.1 - de novo long-read assembly (Flye, Oxford Nanopore)
# In : nanopore_reads.fastq.gz
# Out: FoxMRL8996_AssemblyScaffolds.fasta (draft contigs)
set -euo pipefail

READS="nanopore_reads.fastq.gz"          # [CONFIRM path]
OUTDIR="flye_assembly"
THREADS=32                               # [FILL]

# --nano-hq for high-accuracy ONT (Q20+ chemistry); --nano-raw for older reads. [CONFIRM]
flye --nano-hq "$READS" --out-dir "$OUTDIR" --threads "$THREADS"
# Older Flye (<2.9) needs a genome-size estimate: add  --genome-size 52m   [only if required]

# (optional) ONT consensus polishing - [CONFIRM if used]
# medaka_consensus -i "$READS" -d "$OUTDIR/assembly.fasta" -o medaka -t "$THREADS" -m [FILL model]

cp "$OUTDIR/assembly.fasta" FoxMRL8996_AssemblyScaffolds.fasta
echo "Draft assembly -> FoxMRL8996_AssemblyScaffolds.fasta"
