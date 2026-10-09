# Ship the MEASURED reliability curve, because the one-parameter closed form cannot reproduce it.
# Per-individual (a, w, iota) from the genome-wide held-out-founder simulation, binned on q_hat.
suppressPackageStartupMessages(library(data.table))
L <- "AM52k/lai_sim"
S <- rbindlist(lapply(Sys.glob(file.path(L,"data_gw/chr*/score_rfmix.tsv")), fread))
SS <- c("n","sC","sT","sCT","sCC","sTT"); S[, (SS) := lapply(.SD, as.numeric), .SDcols=SS]
G <- S[, .(q=q[1], qhat=sum(sC)/sum(2*n), n=sum(n), sC=sum(sC), sT=sum(sT),
           sCT=sum(sCT), sCC=sum(sCC), sTT=sum(sTT)), by=.(sample, focal)]
G[, `:=`(a=(sCT-sC*sT/n)/(n-1), w=(sCC-sC^2/n)/(n-1), iota=(sTT-sT^2/n)/(n-1))]
G <- G[is.finite(a)&is.finite(w)&is.finite(iota)]
br <- c(-Inf, 10^seq(-3.6, -0.3, length.out = 28), Inf)
CUR <- G[, .(q_hat=median(qhat), a=mean(a), w=mean(w), v=mean(iota), n=.N),
         by=.(focal, b=cut(qhat, br))][n >= 15][, b := NULL][order(focal, q_hat)]
CUR[, `:=`(rho = a/w, eff_i = a^2/w/pmax(v, 1e-12))]
adm_curves <- as.data.frame(CUR)
attr(adm_curves, "source") <- paste(
  "RFMix 2 against the 5-way HGDP+1kGP reference used in UK Biobank (AFR/EAS/EUR/NAM/SAS);",
  "16,000 held-out-founder mosaics, T = 8 generations, q drawn from the real 52,639 distribution.")
attr(adm_curves, "eps0") <- c(AFR = 1.4e-3, SAS = 4.4e-3)
save(adm_curves, file = "ADMIT/R/sysdata.rda", compress = "xz")
cat(sprintf("shipped %d bins (%s)\n", nrow(adm_curves),
            paste(sprintf("%s %d", names(table(adm_curves$focal)), table(adm_curves$focal)), collapse=", ")))
print(head(adm_curves[adm_curves$focal=="AFR", c("q_hat","rho","n")], 8), digits=3)
