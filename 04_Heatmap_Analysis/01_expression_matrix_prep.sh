#!/usr/bin/env bash
# Stage 04.1 - build the similarity matrix the heatmap consumes.
# NAMING NOTE: despite the filename, this produces a TM-align STRUCTURAL matrix
#   (tm_align_results_all_vs_all.tsv: Protein_A, Protein_B, ..., Max_TM_Score),
#   which is exactly what 02_plot_heatmap.R reads. It is NOT a transcript-
#   expression matrix. If you also have RNA-seq expression, keep it separate. [CONFIRM]
set -euo pipefail

PDB_DIR="PDB_clean"; OUT="tm_align_results_all_vs_all.tsv"
printf 'Protein_A\tProtein_B\tTM1\tTM2\tMax_TM_Score\n' > "$OUT"

pdbs=("$PDB_DIR"/*.pdb); n=${#pdbs[@]}
for ((a=0; a<n; a++)); do
  for ((b=a+1; b<n; b++)); do
    A=$(basename "${pdbs[a]}" .pdb); B=$(basename "${pdbs[b]}" .pdb)
    res=$(TMalign "${pdbs[a]}" "${pdbs[b]}" | grep -E '^TM-score=')
    tm1=$(printf '%s\n' "$res" | sed -n '1s/.*= *\([0-9.]*\).*/\1/p')
    tm2=$(printf '%s\n' "$res" | sed -n '2s/.*= *\([0-9.]*\).*/\1/p')
    max=$(python3 -c "print(max($tm1, $tm2))")
    printf '%s\t%s\t%s\t%s\t%s\n' "$A" "$B" "$tm1" "$tm2" "$max" >> "$OUT"
  done
done
echo "TM-align matrix -> $OUT"
# RIGOR: Max_TM_Score normalizes by the SHORTER chain -> inflates similarity for
#   size-mismatched pairs. State this in Methods; don't over-read cross-family
#   edges that join very different-length proteins.
