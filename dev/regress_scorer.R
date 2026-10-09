# Regression test for adm_score: does it reproduce the SHIPPED curve from the raw .msp?
# The shipped curve was built by a different code path (build_curve.R, from precomputed sufficient
# statistics in score_rfmix.tsv). adm_score goes from the msp and the truth directly. One chromosome
# of the genome-wide arm is enough to compare reliability at matched q_hat.
for (f in list.files("ADMIT/R", pattern = "[.]R$", full.names = TRUE)) source(f)
load("ADMIT/R/sysdata.rda")
suppressPackageStartupMessages(library(data.table))
D  <- "AM52k/lai_sim/data_gw/chr22"
tr <- readRDS(file.path(D, "truth.rds"))
fo <- tr$design$focal == "AFR"
sim <- structure(list(design = as.data.frame(tr$design)[fo, ], site = as.data.frame(tr$site),
                      anc1 = tr$anc1[, fo], anc2 = tr$anc2[, fo],
                      focal = "AFR", baseline = "EUR"), class = "admit_sim")
W <- .read_msp(file.path(D, "rfmix.msp.tsv"), sim, NULL)
W$calls <- W$calls[, rep(fo, each = 2), drop = FALSE]
cur <- adm_score(sim, called = W$calls[W$idx, , drop = FALSE], format = "matrix", min_n = 10L)

# the shipped curve, restricted to the same q_hat range
sh <- adm_curve("AFR")
rng <- range(cur$q_hat[cur$q_hat > 0])
cmp <- data.frame(q_hat = cur$q_hat,
                  rho_adm_score = cur$rho,
                  rho_shipped   = approx(sh$q_hat, sh$rho, xout = cur$q_hat, rule = 2)$y)
cmp <- cmp[cur$q_hat > 0, ]
cat("reliability: adm_score (from the msp, chr22) vs the shipped curve (genome-wide)\n")
print(head(cmp, 12), digits = 3)
cat(sprintf("\nmedian |difference| = %.3f over %d bins\n",
            median(abs(cmp$rho_adm_score - cmp$rho_shipped)), nrow(cmp)))
cat(sprintf("Spearman = %.3f\n", cor(cmp$rho_adm_score, cmp$rho_shipped, method = "spearman")))
