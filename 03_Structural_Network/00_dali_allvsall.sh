#!/usr/bin/env bash
# Stage 03.0 - all-vs-all structural comparison (DaliLite v5).
# Produces the 'ordered' Z-score matrix + dali_mapping_dict.tsv consumed by
# 01_network_construction.R and 03_family_representatives.R.
# Native all-vs-all (matrix) mode - NOT snakedali (query-vs-database).
set -euo pipefail

PDB_DIR="PDB_clean"; DAT_DIR="DAT"; mkdir -p "$DAT_DIR"

# 1. import each PDB -> DAT with a short Dali_ID (P001, P002, ...) + build mapping
printf 'Dali_ID\tJGI_Protein_ID\theader\n' > dali_mapping_dict.tsv
i=0
for pdb in "$PDB_DIR"/*.pdb; do
  i=$((i+1)); id=$(printf 'P%03d' "$i"); base=$(basename "$pdb" .pdb)
  import.pl --pdbfile "$pdb" --pdbid "$id" --dat "$DAT_DIR" --clean
  jgi=$(printf '%s' "$base" | sed -n 's/.*FoxMRL8996_\([0-9]\+\).*/\1/p')
  printf '%s\t%s\t%s\n' "$id" "$jgi" "$base" >> dali_mapping_dict.tsv
done

# 2. native all-vs-all matrix
ls "$DAT_DIR" | sed 's/\.dat$//' > ids.list
dali.pl --matrix --query ids.list --dat1 "$DAT_DIR" --clean
# [CHECK] confirm the matrix output filename on your DaliLite build; the R scripts
#   expect a file named 'ordered' (rename if your build calls it differently).
echo "DALI all-vs-all done: 'ordered' + dali_mapping_dict.tsv"
