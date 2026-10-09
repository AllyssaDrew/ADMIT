#' Split a reference panel into held-out founders and a reduced panel
#'
#' Simulated haplotypes must be built from individuals the caller will **not** see in its reference,
#' or it scores them against themselves and the measured error is meaningless. This is the single
#' easiest thing to get wrong in a calibration, and it fails silently -- accuracy simply looks
#' excellent -- so the split is a function rather than an instruction.
#'
#' @param pop Character vector, one entry per reference *individual*, naming their ancestry.
#' @param frac Fraction of each ancestry to hold out.
#' @param ancestries Which ancestries to draw founders from. Default: all.
#' @return List with `founders` and `reference`, both integer indices into `pop`.
#' @examples
#' pop <- rep(c("AFR", "EUR", "SAS"), each = 100)
#' str(adm_holdout(pop, frac = 0.25))
#' @export
adm_holdout <- function(pop, frac = 0.25, ancestries = unique(pop)) {
  stopifnot(frac > 0, frac < 1)
  idx <- seq_along(pop)
  fnd <- unlist(lapply(ancestries, function(a) {
    i <- idx[pop == a]
    if (length(i) < 8) stop("ancestry '", a, "' has only ", length(i), " individuals")
    sample(i, max(2L, floor(frac * length(i))))
  }), use.names = FALSE)
  list(founders = sort(fnd), reference = sort(setdiff(idx, fnd)))
}

#' Simulate admixed mosaics from held-out founders
#'
#' Each simulated haplotype is a Markov chain along the genetic map: switch points are Poisson with
#' rate `T` per Morgan, each resulting tract takes focal ancestry with probability `q` and the
#' baseline otherwise, and copies a random founder haplotype of that ancestry. Both haplotypes of an
#' individual are drawn independently, so the simulated cohort has no correlation between an
#' individual's two copies.
#'
#' `q` is drawn from `q_dist`, which should be your own cohort's global ancestry. This matters more
#' than it looks: some callers estimate the number of generations from the whole cohort, so a
#' mostly-unadmixed cohort changes how they smooth the tracts of the few admixed individuals. A
#' uniform grid of `q` would hide that.
#'
#' @param hap Integer matrix of alleles, sites in rows, founder *haplotypes* in columns.
#' @param hap_pop Ancestry of each column of `hap`.
#' @param site Data frame of sites with `CHROM`, `POS`, and optionally `ID`, `REF`, `ALT`.
#' @param map Data frame with `pos` and `cM`, the genetic map for this chromosome.
#' @param q_dist Numeric vector to draw `q` from -- your cohort's global focal ancestry.
#' @param focal,baseline Ancestry labels, matching values in `hap_pop`.
#' @param n Number of individuals to simulate.
#' @param T Generations since admixture (Poisson rate per Morgan).
#' @param out_dir If given, writes `sim_query.vcf` there (gzipped and indexed when `bcftools` is on
#'   the path), which is what you then run your own caller on.
#' @return An object of class `admit_sim` carrying the design and the per-site truth.
#' @export
adm_simulate <- function(hap, hap_pop, site, map, q_dist, focal, baseline = "EUR",
                         n = 8000L, T = 8L, out_dir = NULL) {
  stopifnot(nrow(hap) == nrow(site), ncol(hap) == length(hap_pop))
  if (!all(c(focal, baseline) %in% hap_pop))
    stop("`hap_pop` must contain both '", focal, "' and '", baseline, "'")
  M <- approx(map$pos, map$cM, xout = site$POS, rule = 2)$y / 100     # Morgans
  pool <- split(seq_along(hap_pop), hap_pop)
  rng <- range(M); L <- diff(rng)

  one <- function(q) {
    br  <- sort(runif(rpois(1L, T * L), rng[1], rng[2]))
    seg <- findInterval(M, br) + 1L
    anc <- ifelse(runif(length(br) + 1L) < q, focal, baseline)
    src <- vapply(anc, function(a) pool[[a]][sample.int(length(pool[[a]]), 1L)], 1L)
    list(allele = hap[cbind(seq_len(nrow(hap)), src[seg])], focal = anc[seg] == focal)
  }
  q <- sample(q_dist[is.finite(q_dist)], n, replace = TRUE)
  A1 <- A2 <- matrix(0L, nrow(hap), n); F1 <- F2 <- matrix(FALSE, nrow(hap), n)
  for (i in seq_len(n)) {
    a <- one(q[i]); b <- one(q[i])
    A1[, i] <- a$allele; A2[, i] <- b$allele; F1[, i] <- a$focal; F2[, i] <- b$focal
  }
  design <- data.frame(sample = sprintf("sim_%s_%05d", focal, seq_len(n)),
                       focal = focal, q = q, T = T, stringsAsFactors = FALSE)
  out <- structure(list(design = design, site = site, M = M, anc1 = F1, anc2 = F2,
                        focal = focal, baseline = baseline), class = "admit_sim")
  if (!is.null(out_dir)) { dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    .write_vcf(A1, A2, site, design$sample, file.path(out_dir, "sim_query.vcf"))
    saveRDS(out, file.path(out_dir, "truth.rds")) }
  out
}

.write_vcf <- function(A1, A2, site, samples, vf) {
  hdr <- c("##fileformat=VCFv4.2",
           sprintf("##contig=<ID=%s>", unique(site$CHROM)),
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", samples),
                 collapse = "\t"))
  writeLines(hdr, vf)
  body <- cbind(site$CHROM, site$POS,
                if (!is.null(site$ID)) site$ID else ".",
                if (!is.null(site$REF)) site$REF else "A",
                if (!is.null(site$ALT)) site$ALT else "G",
                ".", "PASS", ".", "GT", matrix(paste0(A1, "|", A2), nrow(A1)))
  utils::write.table(body, vf, sep = "\t", row.names = FALSE, col.names = FALSE,
                     quote = FALSE, append = TRUE)
  if (nzchar(Sys.which("bcftools")))
    system2("bcftools", c("view","-Oz","-o", paste0(vf, ".gz"), vf)) &&
      system2("bcftools", c("index","-t", paste0(vf, ".gz")))
  invisible(vf)
}

#' @export
print.admit_sim <- function(x, ...) {
  cat(sprintf("admit_sim: %d individuals, %d sites, focal %s vs %s, T = %d\n",
              nrow(x$design), nrow(x$site), x$focal, x$baseline, x$design$T[1]))
  cat(sprintf("  q drawn from your cohort: median %.4g, %.1f%% above 0.02\n",
              median(x$design$q), 100 * mean(x$design$q > 0.02)))
  invisible(x)
}
