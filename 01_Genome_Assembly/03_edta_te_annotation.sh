#!/usr/bin/env bash
# Stage 01.3 - transposable element annotation (EDTA v2.2.x) on the final assembly.
# In : MRL8996_chromosome_level.fasta (+ CDS to protect genes from masking)
# Out: <genome>.mod.EDTA.TElib.fa      non-redundant TE library
#      <genome>.mod.EDTA.TEanno.gff3   whole-genome annotation (needs --anno 1)
#      <genome>.mod.EDTA.TEanno.sum    per-superfamily summary (%% of genome)
#      <genome>.mod.EDTA.intact.gff3   structurally intact elements
#      <genome>.mod.MAKER.masked       masked genome
set -euo pipefail

GENOME="MRL8996_chromosome_level.fasta"    # [CONFIRM path]
CDS="[FILL FoxMRL8996_CDS.fasta]"          # predicted CDS - strongly recommended (see note)
THREADS=32

# fungal genome -> --species others  (NEVER Rice/Maize: those seed a plant TE library
#   and bias the annotation of a Fusarium genome).
# --anno 1     required for the whole-genome TEanno.gff3 (not just the library).
# --sensitive 1 runs RepeatModeler on the remainder (slower, more complete).
# --cds        lets EDTA exclude real genes -> avoids masking effectors as TEs.
EDTA.pl \
  --genome "$GENOME" \
  --cds "$CDS" \
  --species others \
  --step all \
  --sensitive 1 \
  --anno 1 \
  --evaluate 1 \
  --threads "$THREADS"
# Requires bedtools on PATH, or the density/divergence plots fail silently.

# NOTE: EDTA's default neutral rate (--u 1.3e-8) is rice-derived. If you report LTR
#   insertion AGES, set a Fusarium-appropriate rate via --u, or the ages are mis-dated.

# ---- downstream: effector <-> TE co-localisation (Fusarium LS compartments) ----
# F. oxysporum effectors concentrate in TE-rich lineage-specific regions, so the
# natural link to 02_Effector_Prediction is a coordinate intersect:
# bedtools intersect -a effectors.bed \
#   -b "${GENOME}.mod.EDTA.TEanno.gff3" -wa -wb > effectors_in_TE.tsv
echo "EDTA done -> ${GENOME}.mod.EDTA.TEanno.gff3 / .TElib.fa / .sum"
