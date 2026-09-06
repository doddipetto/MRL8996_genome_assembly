#!/bin/bash
# Fold the Fo47 final effectorome with the local ColabFold install.
# Mirrors the MRL8996 run: colabfold_batch, 3 models, rank_001 retained.
# Resumable: effectors with an existing rank_001 model are skipped.
set -u

FASTA="/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026/data/Fo47_effectors_for_colabfold.fasta"
OUT_DIR="/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026/results/colabfold_Fo47"
COLABFOLD_BIN="/home/usuario/localcolabfold/.pixi/envs/default/bin/colabfold_batch"
PAUSE_SEC=60          # GPU cool-down between sequences, as in the MRL8996 run

mkdir -p "$OUT_DIR" "$OUT_DIR/_split"

# one FASTA per sequence so the run can be resumed per effector
awk '/^>/ { acc = substr($1, 2); f = "'"$OUT_DIR"'/_split/" acc ".fasta" }
     { print > f }' "$FASTA"

for fa in "$OUT_DIR"/_split/*.fasta; do
    base=$(basename "$fa" .fasta)
    if ls "$OUT_DIR"/${base}*rank_001* 1> /dev/null 2>&1; then
        echo "[$base] model already present, skipping"
        continue
    fi
    echo "[$base] folding (3 models)"
    "$COLABFOLD_BIN" --num-models 3 "$fa" "$OUT_DIR"
    sleep $PAUSE_SEC
done

echo "Fo47 effectorome folding finished; models in $OUT_DIR"
