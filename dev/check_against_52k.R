# Does the closed form reproduce what the truth simulation measured on the same cohort?
# This is the check that decides whether tier 2 (no simulation) is worth shipping.
for (f in list.files("ADMIT/R", full.names = TRUE)) source(f)
cv <- read.delim("AM52k/data/cohort_covariates_52k.tsv")
cat("cohort:", nrow(cv), "individuals\n\n")

for (focal in c("AFR","SAS")) {
  Q <- cv[, c("AFR","EAS","EUR","NAM","SAS")]
  # genome-wide mean focal local dosage is not stored per individual; the global Q IS that mean,
  # so for this check we use it as the dosage (they are the same quantity by construction)
  e <- adm_eps0(Q, local_dosage = 2 * Q[[focal]], focal = focal, pure_cut = 0.99)
  print(e); cat("\n")
  th <- adm_threshold(cv[[focal]], eps0 = e)
  print(th); cat("\n--------------------------------------------------\n")
}
