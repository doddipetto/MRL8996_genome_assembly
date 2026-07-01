#!/usr/bin/env bash
# Stage 02.1 - secretome prediction (SignalP 6) + MATURE-sequence extraction.
# In : predicted proteome (FoxMRL8996_proteins.fasta)   [CONFIRM name]
# Out: FoxMRL8996_secreted_mature.fasta (signal peptide removed)
set -euo pipefail

PROTEOME="[FILL FoxMRL8996_proteins.fasta]"
OUTDIR="signalp6_out"

signalp6 --fastafile "$PROTEOME" --organism eukarya \
         --output_dir "$OUTDIR" --format none --mode fast
# SignalP6 writes processed_entries.fasta = MATURE sequences (SP cleaved).
# CRITICAL downstream: AMAPEC (02.5) must fold MATURE seqs, so keep these.
cp "$OUTDIR/processed_entries.fasta" FoxMRL8996_secreted_mature.fasta
echo "Secreted mature proteins -> FoxMRL8996_secreted_mature.fasta"
