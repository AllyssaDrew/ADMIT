for (f in list.files("ADMIT/R", full.names = TRUE)) source(f)
cv <- read.delim("AM52k/data/cohort_covariates_52k.tsv")
cat("闭式解对 eps0 的敏感度（焦点 AFR，实测最优 = 4.4e-3，实测效率 = 0.784）\n\n")
cat(sprintf("%-12s %-12s %-10s %-10s\n", "eps0", "detection", "eff", "keeps"))
for (e in c(5.87e-4, 1.0e-3, 1.4e-3, 2e-3, 3e-3, 5e-3)) {
  th <- adm_threshold(cv$AFR, eps0 = e)
  cat(sprintf("%-12.3g %-12.3g %-10.3f %-10s\n", e, th$detection$cut, th$detection$eff,
              format(th$detection$n, big.mark=",")))
}
cat("\n闭式信度 vs 真值模拟实测（RESULTS 3: 0.14 @ 0.001, 0.74 @ 0.007, 0.95 @ >0.05）\n")
for (e in c(5.87e-4, 1.4e-3)) {
  cat(sprintf("  eps0=%.3g : r(0.001)=%.2f  r(0.007)=%.2f  r(0.05)=%.2f\n", e,
      adm_reliability(max((0.001-e),0), e), adm_reliability(max((0.007-e),0), e),
      adm_reliability(max((0.05-e),0), e)))
}
