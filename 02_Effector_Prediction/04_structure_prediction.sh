#!/usr/bin/env bash
# Stage 02.4 - 3D structure prediction of the FINAL effectorome.
# CRITICAL: fold MATURE sequences (SP already removed). Folding WITH the signal
#   peptide distorts geometry and biases AMAPEC (per AMAPEC docs).
# Out: one PDB per effector in PDB_clean/  (shared input: AMAPEC, TM-align, DALI)
set -euo pipefail

EFF_FASTA="FoxMRL8996_Final_Effectorome.fasta"   # from 03 (mature, non-TM)
RAW="colabfold_output"; mkdir -p "$RAW" PDB_clean

# Option A - ColabFold (better MSA; GPU). RTX 4070 SUPER 12 GB handles effector-sized proteins.
colabfold_batch "$EFF_FASTA" "$RAW" --num-recycle 3
# keep rank_001 per protein, strip ColabFold suffix so IDs match canon() in the heatmap:
for p in "$RAW"/*_rank_001_*.pdb; do
  [ -e "$p" ] || continue
  base=$(basename "$p" | sed -E 's/_(un)?relaxed_rank.*//')
  cp "$p" "PDB_clean/${base}.pdb"
done

# Option B - ESMFold (faster, lower compute, minor accuracy drop) - [CONFIRM which you used]
# python esmfold_batch.py -i "$EFF_FASTA" -o PDB_clean

echo "Per-effector PDBs -> PDB_clean/"
