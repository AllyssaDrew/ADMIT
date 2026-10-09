# Tier 3 scored on data the package never saw: the gnomAD-panel chr22 arm, whose truth.rds and
# RFMix .msp came from the original pipeline rather than from adm_simulate().
for (f in list.files("ADMIT/R", pattern = "[.]R$", full.names = TRUE)) source(f)
load("ADMIT/R/sysdata.rda")
D  <- "AM52k/lai_sim/data_gnomad"
tr <- readRDS(file.path(D, "truth.rds"))
fo <- tr$design$focal == "AFR"
sim <- structure(list(design = as.data.frame(tr$design)[fo, ], site = as.data.frame(tr$site),
                      anc1 = tr$anc1[, fo], anc2 = tr$anc2[, fo],
                      focal = "AFR", baseline = "EUR"), class = "admit_sim")
cat(sprintf("%d focal-AFR individuals, %d sites\n", nrow(sim$design), nrow(sim$site)))

msp <- file.path(D, "out_rfmix_way5/out.msp.tsv")
W <- .read_msp(msp, sim, NULL)
W$calls <- W$calls[, rep(fo, each = 2), drop = FALSE]      # same individuals, two columns each
cur <- adm_score(sim, called = W$calls[W$idx, , drop = FALSE], format = "matrix")

cat("\ncurve measured by tier 3 (chr22, gnomAD 5-way panel):\n")
print(head(cur[, c("q_hat","rho","n")], 7), digits = 3)
cat(sprintf("eps0 implied: %.3g\n", attr(cur, "eps0")))
cat("\nshipped curve (genome-wide, production 5-way panel):\n")
sh <- adm_curve("AFR"); print(head(sh[, c("q_hat","rho","n")], 7), digits = 3)
cv <- read.delim("AM52k/data/cohort_covariates_52k.tsv")
cat(sprintf("\nwindow from tier-3 curve  : %.3g (eff %.3f)\n",
    adm_threshold(cv$AFR, curve = cur)$detection$cut, adm_threshold(cv$AFR, curve = cur)$detection$eff))
cat(sprintf("window from shipped curve : %.3g (eff %.3f)\n",
    adm_threshold(cv$AFR, curve = sh)$detection$cut, adm_threshold(cv$AFR, curve = sh)$detection$eff))
