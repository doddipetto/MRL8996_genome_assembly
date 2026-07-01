#!/usr/bin/env bash
# Stage 02.2 - effector prediction (EffectorP 3.0), intersect with secretome,
#              then transmembrane screen (DeepTMHMM via BioLib).
set -euo pipefail

SECRETED_MATURE="FoxMRL8996_secreted_mature.fasta"   # from 02.1

# EffectorP 3.0 on the secreted set. [CONFIRM invocation for your install]
python EffectorP.py -i "$SECRETED_MATURE" -o effectorp3_out.txt -E effector_candidates.fasta

# Because EffectorP ran on the SignalP-secreted set, its output already IS the
# SignalP ∩ EffectorP intersect:
cp effector_candidates.fasta SignalP_EffectorP_intersect.fasta

# Transmembrane screen (DeepTMHMM, BioLib CLI):
# pip install -U pybiolib
biolib run DTU/DeepTMHMM --fasta SignalP_EffectorP_intersect.fasta
# -> biolib_results/predicted_topologies.3line  (consumed by 03_parse_effectors.py)
echo "Now run 03_parse_effectors.py to drop TM-containing proteins."
