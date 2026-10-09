# Can a reliability curve measured for one caller be rescaled to another by the ratio of their eps0?
#
# The closed form says yes: r(q) = 2q/(2q+eps0) depends on q only through q/eps0, so plotting r
# against q/eps0 should collapse any two callers onto one line. RFMix and FLARE were run on the SAME
# 16,000 mosaics, so the question is decidable rather than arguable.
#
# Three things are measured, in increasing order of what a user cares about:
#   1. eps0 for each caller, at true q = 0
#   2. whether the two measured curves collapse under q -> q/eps0
#   3. the operational cost: the window you pick from FLARE's own curve, versus from the shipped
#      RFMix curve rescaled, versus from the shipped RFMix curve used naively
suppressPackageStartupMessages(library(data.table))
for (f in list.files("ADMIT/R", pattern = "[.]R$", full.names = TRUE)) source(f)
L <- "AM52k/lai_sim"

load_tool <- function(tag) {
  S <- rbindlist(lapply(Sys.glob(file.path(L, sprintf("data_gw/chr*/score_%s.tsv", tag))), fread))
  SS <- c("n","sC","sT","sCT","sCC","sTT"); S[, (SS) := lapply(.SD, as.numeric), .SDcols = SS]
  G <- S[, .(q = q[1], qhat = sum(sC)/sum(2*n), n = sum(n), sC = sum(sC), sT = sum(sT),
             sCT = sum(sCT), sCC = sum(sCC), sTT = sum(sTT)), by = .(sample, focal)]
  G[, `:=`(a = (sCT - sC*sT/n)/(n-1), w = (sCC - sC^2/n)/(n-1), iota = (sTT - sT^2/n)/(n-1))]
  G[focal == "AFR" & is.finite(a) & is.finite(w) & is.finite(iota)][, tool := tag][]
}
R <- load_tool("rfmix"); F <- load_tool("flare")
cat(sprintf("same mosaics: RFMix %d individuals, FLARE %d; shared samples %d\n\n",
            nrow(R), nrow(F), length(intersect(R$sample, F$sample))))

## ---- 1. eps0 at true q = 0 ---------------------------------------------------------------------
e <- c(rfmix = R[q == 0, mean(qhat)], flare = F[q == 0, mean(qhat)])
cat(sprintf("eps0 at true q = 0 :  RFMix %.3g   FLARE %.3g   ratio %.2f\n",
            e["rfmix"], e["flare"], e["rfmix"]/e["flare"]))
cat(sprintf("calling gain lambda :  RFMix %.3f   FLARE %.3f\n\n",
            R[q > .05, sum(a)/sum(iota)], F[q > .05, sum(a)/sum(iota)]))

## ---- 2. do the curves collapse under q -> q/eps0 ? ---------------------------------------------
bin <- function(G, eps) {
  G <- copy(G)[q > 0]
  G[, s := q / eps]
  br <- 10^seq(-1.2, 3.2, length.out = 16)
  G[, b := cut(s, br)]
  G[!is.na(b), .(s = median(s), rho = sum(a)/sum(w), n = .N), by = b][n >= 20][order(s)]
}
BR <- bin(R, e["rfmix"]); BF <- bin(F, e["flare"])
M <- merge(BR[, .(s_r = s, rho_r = rho)], BF[, .(s_f = s, rho_f = rho)],
           by.x = "s_r", by.y = "s_f", all = FALSE)
if (!nrow(M)) {   # bins rarely align exactly; interpolate FLARE onto RFMix's scaled grid
  M <- data.table(s = BR$s, rho_r = BR$rho,
                  rho_f = approx(BF$s, BF$rho, xout = BR$s, rule = 2)$y)
} else setnames(M, c("s","rho_r","rho_f"))
cat("reliability against the SCALED axis q/eps0 -- if rescaling worked these two columns agree:\n")
print(M[, .(q_over_eps0 = round(s, 2), rho_RFMix = round(rho_r, 3),
            rho_FLARE = round(rho_f, 3), diff = round(rho_f - rho_r, 3))])
cat(sprintf("\n  median |difference| across the shared range: %.3f\n", median(abs(M$rho_f - M$rho_r))))

## ---- 3. the operational cost -------------------------------------------------------------------
mkcurve <- function(G) {
  br <- c(-Inf, 10^seq(-3.6, -0.3, length.out = 28), Inf)
  d <- G[, .(q_hat = median(qhat), a = mean(a), w = mean(w), v = mean(iota), n = .N),
         by = .(b = cut(qhat, br))][n >= 15][order(q_hat)]
  as.data.frame(d[, .(q_hat, a, w, v, n)])
}
cv <- read.delim("AM52k/data/cohort_covariates_52k.tsv")
cR <- mkcurve(R); cF <- mkcurve(F)
cRs <- cR; cRs$q_hat <- cR$q_hat * e["flare"] / e["rfmix"]      # rescale RFMix -> FLARE by eps0
arms <- list("FLARE's own curve (truth)" = cF,
             "shipped RFMix curve, rescaled by eps0" = cRs,
             "shipped RFMix curve, used naively"     = cR)
cat("\nwindow a FLARE user would pick, three ways:\n")
for (nm in names(arms)) {
  th <- adm_threshold(cv$AFR, curve = arms[[nm]])
  cat(sprintf("  %-40s detection %-9.3g eff %.3f   keeps %s\n", nm,
              th$detection$cut, th$detection$eff, format(th$detection$n, big.mark = ",")))
}

## ---- 4. what the rescaled choice actually COSTS, scored on FLARE's own curve -------------------
truth <- adm_efficiency(cv$AFR, curve = cF)
best  <- truth[which.max(truth$eff), ]
pick  <- function(cur) adm_threshold(cv$AFR, curve = cur)$detection$cut
at    <- function(c0) truth$eff[which.min(abs(log(truth$cut) - log(c0)))]
cat("\nefficiency ACTUALLY attained, every cut scored on FLARE's own curve:\n")
cat(sprintf("  best available                          %.3f  at %.3g\n", best$eff, best$cut))
for (nm in names(arms)[-1])
  cat(sprintf("  cut chosen from %-24s %.3f  at %.3g   (costs %.3f)\n",
              sub("shipped RFMix curve, ", "", nm), at(pick(arms[[nm]])), pick(arms[[nm]]),
              best$eff - at(pick(arms[[nm]]))))
