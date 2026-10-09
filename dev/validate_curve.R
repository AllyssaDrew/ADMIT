for (f in list.files("ADMIT/R", pattern = "[.]R$", full.names = TRUE)) source(f)
load("ADMIT/R/sysdata.rda")
cv <- read.delim("AM52k/data/cohort_covariates_52k.tsv")
cat("measured-curve mode vs what we measured directly on this cohort\n")
cat("  AFR measured: optimum 4.4e-3, efficiency 0.784 | SAS measured: 2.6e-3, 0.477\n\n")
for (fo in c("AFR", "SAS")) {
  th <- adm_threshold(cv[[fo]], curve = fo)
  cat(sprintf("  %-4s detection %-9.3g eff %.3f | estimation %-9.3g eff %.3f | keeps %s\n",
      fo, th$detection$cut, th$detection$eff, th$estimation$cut, th$estimation$eff,
      format(th$detection$n, big.mark = ",")))
}
