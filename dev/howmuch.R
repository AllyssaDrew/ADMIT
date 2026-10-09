# 为什么 chr22 不够：把"每个 bin 里有多少条真实的焦点片段"算出来
suppressPackageStartupMessages(library(data.table))
for (f in list.files("ADMIT/R", pattern="[.]R$", full.names=TRUE)) source(f)
load("ADMIT/R/sysdata.rda")

# 各染色体的遗传长度（Morgan），从 RFMix 用的同一张图
GM <- "/proj/xihaoli_lab/users/quansun/local_ans_pipeline/rfmix_files/genetic_map_hg38_withX.sorted.txt"
gm <- fread(GM, header=TRUE, showProgress=FALSE); setnames(gm, 1:3, c("chr","pos","cM")); gm[, cM := as.numeric(cM)]
gm <- gm[chr %in% paste0("chr", 1:22)]
Lchr <- gm[, .(M = (max(cM)-min(cM))/100), by=chr]
L_gw <- sum(Lchr$M); L_22 <- Lchr[chr=="chr22", M]
cat(sprintf("基因组遗传长度 %.1f Morgan；chr22 %.2f Morgan（%.1f%%，即 %.0f 倍差距）\n\n",
            L_gw, L_22, 100*L_22/L_gw, L_gw/L_22))

# 一个个体在长度 L 上期望有多少条【焦点】片段 = T * L * q
cat("单个个体期望的焦点片段数  (T = 8)\n")
cat(sprintf("%-10s %-14s %-14s\n", "q", "全基因组", "chr22"))
for (q in c(5e-4, 1e-3, 5e-3, 2e-2, 5e-2)) {
  cat(sprintf("%-10.4g %-14.2f %-14.3f\n", q, 8*L_gw*q, 8*L_22*q))
}
# 一个 bin 里的总片段数 = n_bin * T * L * q  —— 这才是决定 a、w 均值稳不稳的量
sh <- adm_curve("AFR")
cat("\nship 的曲线里每个 bin 的总焦点片段数（全基因组 vs 若只用 chr22）\n")
cat(sprintf("%-10s %-7s %-14s %-14s %-8s\n", "q_hat", "n_bin", "全基因组", "chr22", "信度"))
for (i in seq_len(min(10, nrow(sh)))) {
  q <- sh$q_hat[i]; n <- sh$n[i]
  cat(sprintf("%-10.4g %-7d %-14.1f %-14.2f %-8.3f\n", q, n, n*8*L_gw*q, n*8*L_22*q, sh$rho[i]))
}
cat("\n经验法则：bin 内总片段数 ≳ 30 才稳；chr22 上要到 q_hat ~ 0.025 才够，正是实测一致性断开的位置。\n")
