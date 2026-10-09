#!/bin/bash
set -uo pipefail
module purge >/dev/null 2>&1; module load intel/2025.0.1 >/dev/null 2>&1
source "$MKLROOT/env/vars.sh" intel64 >/dev/null 2>&1
export R_LIBS_USER=/proj/xihaoli_lab/users/xihaoli/R-4.5.0-MKL
F=/proj/xihaoli_lab/users/ganzeyu/UKB-admixture_mapping/local_ans_pipeline/Admixture_Mapping/Framework
LIB=$(mktemp -d)
$HOME/R-4.5.0/bin/R CMD INSTALL -l "$LIB" "$F/ADMIT" >/dev/null 2>&1
cd "$F/ADMIT"
R_LIBS="$LIB:$R_LIBS_USER" $HOME/R-4.5.0/bin/Rscript -e '
  library(testthat); library(ADMIT); test_dir("tests/testthat", reporter="summary")' 2>&1 | tail -40
