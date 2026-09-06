#!/usr/bin/env python3
"""
06_prepare_Fo47_folding.py
==========================
Prepares the Fo47 arm of the cross-strain structural comparison. Fol4287 has
AlphaFold2 models already (Fol4287_genome_paper/colabfold_results), Fo47 does
not, so its 393 final effectors have to be folded locally before script 05 can
be re-run against them.

This script only WRITES the inputs - it does not fold anything:
    * Fo47_effectors_for_colabfold.fasta   sanitised headers (accession only)
    * run_colabfold_Fo47.sh                launcher for the local GPU
    * Fo47_header_map.tsv                  accession -> original FASTA header

The folding itself is the heavy step and belongs on the RTX 4070 SUPER:

    bash Revision_2026/scripts/run_colabfold_Fo47.sh

with the same settings used for MRL8996 (colabfold_batch, --num-models 3,
rank_001 kept). Expect roughly 1-3 min per effector for short secreted
proteins, i.e. ~10-18 h for 393 sequences, and ~25 GB of output. The launcher
skips any effector whose rank_001 model already exists, so it can be stopped
and resumed.

Once models are in Revision_2026/results/colabfold_Fo47/, re-run script 05
with --targets fo47 to add the Fo47 comparison (the analysis code is identical;
only the target directory and the accession regex change).

Usage: python3 06_prepare_Fo47_folding.py
"""

from pathlib import Path

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
FO47_FASTA = (ADODDI / "F012_HiC_data" / "final_results" / "Effectorome_Final_Results"
              / "Fo47_final_effectors.fasta")
COLABFOLD_BIN = "/home/usuario/localcolabfold/.pixi/envs/default/bin/colabfold_batch"

OUT_FASTA = REV / "data" / "Fo47_effectors_for_colabfold.fasta"
MAP_OUT = REV / "data" / "Fo47_header_map.tsv"
RUNNER = REV / "scripts" / "run_colabfold_Fo47.sh"
MODEL_DIR = REV / "results" / "colabfold_Fo47"

RUNNER_TEMPLATE = """#!/bin/bash
# Fold the Fo47 final effectorome with the local ColabFold install.
# Mirrors the MRL8996 run: colabfold_batch, 3 models, rank_001 retained.
# Resumable: effectors with an existing rank_001 model are skipped.
set -u

FASTA="{fasta}"
OUT_DIR="{out}"
COLABFOLD_BIN="{bin}"
PAUSE_SEC=60          # GPU cool-down between sequences, as in the MRL8996 run

mkdir -p "$OUT_DIR" "$OUT_DIR/_split"

# one FASTA per sequence so the run can be resumed per effector
awk '/^>/ {{ acc = substr($1, 2); f = "'"$OUT_DIR"'/_split/" acc ".fasta" }}
     {{ print > f }}' "$FASTA"

for fa in "$OUT_DIR"/_split/*.fasta; do
    base=$(basename "$fa" .fasta)
    if ls "$OUT_DIR"/${{base}}*rank_001* 1> /dev/null 2>&1; then
        echo "[$base] model already present, skipping"
        continue
    fi
    echo "[$base] folding (3 models)"
    "$COLABFOLD_BIN" --num-models 3 "$fa" "$OUT_DIR"
    sleep $PAUSE_SEC
done

echo "Fo47 effectorome folding finished; models in $OUT_DIR"
"""


def main():
    (REV / "data").mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    n = 0
    with open(FO47_FASTA) as fh, open(OUT_FASTA, "w") as out, open(MAP_OUT, "w") as mp:
        mp.write("accession\toriginal_header\n")
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                acc = line[1:].split()[0]
                mp.write("%s\t%s\n" % (acc, line[1:]))
                out.write(">%s\n" % acc)
                n += 1
            elif line:
                out.write(line + "\n")

    RUNNER.write_text(RUNNER_TEMPLATE.format(
        fasta=OUT_FASTA, out=MODEL_DIR, bin=COLABFOLD_BIN))
    RUNNER.chmod(0o755)

    print("sequences written: %d -> %s" % (n, OUT_FASTA))
    print("launcher: %s" % RUNNER)
    print("models will be written to: %s" % MODEL_DIR)
    if not Path(COLABFOLD_BIN).exists():
        print("NOTE: colabfold_batch not found at %s - edit COLABFOLD_BIN in the "
              "launcher before running." % COLABFOLD_BIN)


if __name__ == "__main__":
    main()
