#!/bin/bash
RSCRIPT="C:/Program Files/R/R-4.5.2/bin/Rscript.exe"
WD="C:/Users/stat-user/iCloudDrive/Desktop/EBMR"

echo "=== Step 1: Installing EBMRalgorithmFast4 ==="
"$RSCRIPT" -e "setwd('$WD'); devtools::install('EBMRalgorithmFast4', quiet=TRUE, upgrade='never')"
echo "Exit code: $?"

echo "=== Step 2: Running compare_se1_se2.R ==="
"$RSCRIPT" "$WD/Simulation_Studies/compare_se1_se2.R"
echo "Exit code: $?"
