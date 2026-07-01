#!/usr/bin/env bash
# Stage 02.5 - antimicrobial-effector prediction (AMAPEC).
# In : PDB_clean/ (structures from MATURE sequences - see 02.4)
# Out: prediction.csv (Protein ID, Prediction, Probability, pLDDT)
set -euo pipefail

# git clone https://github.com/fantin-mesny/amapec && chmod +x amapec/amapec
# conda env create -f amapec/environment.yml
conda activate amapec
amapec -i PDB_clean -o AMAPEC_Effectorome -t 40

# [CHECK] AMAPEC v1.0 writes  AMAPEC_Effectorome/prediction.csv
#   Your heatmap (04) reads   AMAPEC_Effectorome/results/prediction.csv
#   Make the two agree: either move the file or edit amapec_file in 02_plot_heatmap.R.
echo "AMAPEC -> AMAPEC_Effectorome/prediction.csv"
