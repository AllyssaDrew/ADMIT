#!/bin/bash
# R CMD build + check for ADMIT, on a compute node, with the MKL environment rr sets up.
# Vignettes ARE built here: a package with a vignette source and no inst/doc is a check WARNING.
set -uo pipefail
module purge >/dev/null 2>&1; module load intel/2025.0.1 >/dev/null 2>&1
source "$MKLROOT/env/vars.sh" intel64 >/dev/null 2>&1
export R_LIBS_USER=/proj/xihaoli_lab/users/xihaoli/R-4.5.0-MKL
R=$HOME/R-4.5.0/bin/R
PKG=/proj/xihaoli_lab/users/ganzeyu/UKB-admixture_mapping/local_ans_pipeline/Admixture_Mapping/Framework/ADMIT
WORK=$(mktemp -d); cd "$WORK"
"$R" CMD build "$PKG" 2>&1 | tail -3
TGZ=$(ls ADMIT_*.tar.gz 2>/dev/null | head -1)
[ -z "$TGZ" ] && { echo "BUILD FAILED"; exit 1; }
"$R" CMD check --no-manual "$TGZ" >/dev/null 2>&1
echo "=================================================="
grep -E 'NOTE|WARNING|ERROR|^Status' ADMIT.Rcheck/00check.log | sed 's/^/  /'
echo "--- 细节 ---"
awk '/\.\.\. (NOTE|WARNING|ERROR)/{p=3} p>0{print; p--}' ADMIT.Rcheck/00check.log | head -30
