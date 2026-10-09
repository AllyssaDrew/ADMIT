lib <- tempfile(); dir.create(lib)
try(install.packages("ADMIT", repos = NULL, type = "source", lib = lib,
                     INSTALL_opts = "--no-docs"), silent = TRUE)
ok <- "ADMIT" %in% rownames(installed.packages(lib.loc = lib))
cat("install:", if (ok) "OK" else "FAILED", "\n")
if (ok) {
  library(ADMIT, lib.loc = lib)
  cat("exported:", paste(ls("package:ADMIT"), collapse = ", "), "\n")
  cv <- read.delim("AM52k/data/cohort_covariates_52k.tsv")
  Q <- cv[, c("AFR","EAS","EUR","NAM","SAS")]
  e  <- adm_eps0(Q, 2*Q$AFR, focal = "AFR")
  th <- adm_threshold(cv$AFR, curve = "AFR", eps0 = e)
  cat("end-to-end:", sprintf("eps0 %.3g -> cut %.3g, eff %.3f", e$eps0, th$detection$cut,
      th$detection$eff), "\n")
  cat("adm_rescale:", class(try(adm_rescale(adm_curve("AFR"), 1.33e-3, 3.19e-3), silent=TRUE))[1], "\n")
  cat("plot():", class(try(pdf(NULL) -> z, silent=TRUE))[1], "\n")
}
cat("man/ pages:", length(list.files("ADMIT/man")), "\n")
cat("tests/testthat:", length(list.files("ADMIT/tests/testthat")), "\n")
