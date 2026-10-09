#!/bin/bash
# R CMD build + check, on a compute node, with the MKL environment this cluster needs.
set -uo pipefail
module purge >/dev/null 2>&1; module load intel/2025.0.1 >/dev/null 2>&1
source "$MKLROOT/env/vars.sh" intel64 >/dev/null 2>&1
export R_LIBS_USER=/proj/xihaoli_lab/users/xihaoli/R-4.5.0-MKL
R=$HOME/R-4.5.0/bin/R
PKG="$(cd "$(dirname "$0")/.." && pwd)"
WORK=$(mktemp -d); cd "$WORK"
"$R" CMD build "$PKG" 2>&1 | tail -2
TGZ=$(ls ADMIT_*.tar.gz 2>/dev/null | head -1)
[ -z "$TGZ" ] && { echo "BUILD FAILED"; exit 1; }
"$R" CMD check --no-manual "$TGZ" >/dev/null 2>&1
grep -E 'NOTE|WARNING|ERROR|^Status' ADMIT.Rcheck/00check.log | sed 's/^/  /'
