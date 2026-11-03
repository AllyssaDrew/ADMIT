# Build the shipped curves: the original RFMix pair plus the FLARE curves measured in 2026-11.
# Run from the ADMIT package root with the simulation outputs available.
suppressPackageStartupMessages(library(data.table))
S <- "/proj/xihaoli_lab/users/ganzeyu/UKB-admixture_mapping/local_ans_pipeline/Admixture_Mapping/Framework/AM52k/lai_sim"
AN5 <- c("AFR", "EUR", "EAS", "AMR", "SAS")

## ---- the RFMix pair, unchanged (holdout-founder simulation, data_gw) ----
old <- rbindlist(lapply(Sys.glob(file.path(S, "data_gw/chr*/score_rfmix.tsv")), fread))
SS <- c("n", "sC", "sT", "sCT", "sCC", "sTT"); old[, (SS) := lapply(.SD, as.numeric), .SDcols = SS]
Go <- old[, .(q = q[1], qhat = sum(sC) / sum(2 * n), n = sum(n), sC = sum(sC), sT = sum(sT),
              sCT = sum(sCT), sCC = sum(sCC), sTT = sum(sTT)), by = .(sample, focal)]
Go[, `:=`(a = (sCT - sC * sT / n) / (n - 1), w = (sCC - sC^2 / n) / (n - 1), v = (sTT - sT^2 / n) / (n - 1))]
br <- c(-Inf, 10^seq(-3.6, -0.3, length.out = 28), Inf)
RF <- Go[is.finite(a) & is.finite(w) & is.finite(v), .(q_hat = median(qhat), a = mean(a), w = mean(w), v = mean(v), n = .N),
         by = .(focal, b = cut(qhat, br))][n >= 15][, b := NULL][order(focal, q_hat)]
RF[, `:=`(caller = "RFMix", background = "European", rho = a / w)]
e_rf <- Go[q == 0, .(eps0 = mean(qhat)), by = focal]

## ---- the FLARE curves (AoU-style reference; arms A and B) ----
FL <- fread(file.path(S, "data_aou/final_curves.tsv"))[cohort == "AoU" & focal == "AFR"]
FL[, `:=`(caller = "FLARE", background = "European")]
SB <- fread(file.path(S, "data_aou/final_sas_by_background.tsv"))
SB[, `:=`(focal = "SAS", caller = "FLARE", background = bg)][, bg := NULL]
G <- fread(file.path(S, "data_aou/final_individuals.tsv"))

## ---- the American curves (arm C) ----
## Only the Indigenous-American estimand is shipped. The other one built, AMR_PANEL, takes the whole AMR
## reference class as truth -- mostly admixed Latin American samples carrying about half European ancestry --
## and its reliability FALLS as true ancestry rises (0.83 at q = 0.05 down to 0.40 at 0.39), because a Latino
## founder's European segments are counted as American truth and then correctly called European. That is an
## incoherent estimand, not a measurement of the caller, so it is documented and not shipped.
AM <- fread(file.path(S, "data_aou3/amr_curves.tsv"))[focal == "AMR_NAM"]
AM[, `:=`(focal = "AMR", caller = "FLARE",
          background = fifelse(background == "AFR", "African", "European"))]

## eps0 for each curve, computed exactly as adm_eps0() does on real data: the mean called focal ancestry among
## individuals whose largest called NON-focal ancestry exceeds 0.99. The unselected mean over all q = 0 individuals
## is recorded alongside, because the two differ and only the first is comparable with adm_eps0().
eps_pair <- function(d, f) { oth <- as.matrix(d[, paste0("len_", setdiff(AN5, f)), with = FALSE]) / d$len_tot
  k <- apply(oth, 1, max) > 0.99
  c(selected = if (sum(k) >= 20) mean(d$qhat[k]) else NA_real_, unselected = mean(d$qhat), n_sel = sum(k)) }
E <- rbind(
  data.table(caller = "RFMix", focal = "AFR", background = "European", eps0 = e_rf[focal == "AFR", eps0]),
  data.table(caller = "RFMix", focal = "SAS", background = "European", eps0 = e_rf[focal == "SAS", eps0]),
  rbindlist(lapply(list(c("AFR", "European"), c("SAS", "European"), c("SAS", "African")), function(z) {
    d <- G[focal == z[1] & bgclass == z[2] & q == 0]; p <- eps_pair(d, z[1])
    data.table(caller = "FLARE", focal = z[1], background = z[2],
               eps0 = ifelse(is.na(p["selected"]), p["unselected"], p["selected"]),
               eps0_unselected = p["unselected"], n_ref = p["n_sel"]) })), fill = TRUE)
print(E)

adm_curves <- as.data.frame(rbind(
  RF[, .(focal, caller, background, q_hat, a, w, v, n, rho)],
  FL[, .(focal, caller, background, q_hat, a, w, v, n, rho)],
  SB[, .(focal, caller, background, q_hat, a, w, v, n, rho)],
  AM[, .(focal, caller, background, q_hat, a, w, v, n, rho)]))
eps <- setNames(E$eps0, paste(E$caller, E$focal, E$background, sep = "|"))
## The American curves carry NO eps0, deliberately. adm_eps0's reference rule keeps individuals who are >99% of
## one non-focal ancestry; in this simulation NONE of the truly-pure American-focal individuals qualify, against
## 0.7% in All of Us and 9.1% for the African arm. Two numbers obtained by selecting to incomparable depths
## cannot be divided to give a rescale factor -- doing so would compress this curve 23-fold, which is the same
## error that put the American window at a reliability of zero in the first place. Use these curves as measured.
eps[c("FLARE|AMR|European", "FLARE|AMR|African")] <- NA_real_
attr(adm_curves, "eps0") <- eps
attr(adm_curves, "source") <- paste(
  "RFMix 2 and FLARE 0.6.0 against HGDP+1kGP references, scored on genome-wide mosaics of held-out founders.",
  "RFMix: the 5-way 446-sample reference used in UK Biobank, 3,600 mosaics per focal ancestry.",
  "FLARE: a 3,201-sample 5-way reference (AFR/AMR/EAS/EUR/SAS) matching the All of Us setup, 16,000 mosaics,",
  "founders held out in 5 rotating folds with their 1000 Genomes relatives removed, T = 8, q over the whole 0-1 range.",
  "Validated against the real 52,639-participant UK Biobank cohort run through the same FLARE setup.")
save(adm_curves, file = "R/sysdata.rda", compress = "xz")
cat(sprintf("\nshipped %d rows\n", nrow(adm_curves)))
print(as.data.table(adm_curves)[, .(bins = .N, min_q = signif(min(q_hat), 3), max_rho = round(max(rho), 3)),
                                by = .(caller, focal, background)])
