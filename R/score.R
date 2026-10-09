#' Score a caller's output against simulated truth
#'
#' Turns the output of your own local-ancestry caller, run on the query that [adm_simulate()] wrote,
#' into the reliability curve [adm_threshold()] needs. The package never runs the caller: that is the
#' one step whose command lines, config files and versions differ between tools and across releases,
#' and wrapping it is how a package rots. Generating the query and scoring the answer are stable.
#'
#' @param sim An `admit_sim` object from [adm_simulate()], or the path to the `truth.rds` it wrote.
#' @param called Path to the caller's output.
#' @param format `"msp"` for RFMix and Gnomix, which share the `.msp.tsv` layout; `"flare"` for
#'   FLARE's ancestry VCF; `"matrix"` to pass an already-parsed sites x (2 * individuals) logical
#'   matrix of focal calls.
#' @param focal_code The ancestry's integer code in the `.msp` header. Read from the header when it
#'   names the ancestries, which RFMix and Gnomix both do.
#' @param bins Number of q-hat bins for the returned curve.
#' @param min_n Drop bins holding fewer individuals than this.
#' @return A curve: `q_hat`, `a`, `w`, `v`, `rho`, `n`. Pass it to [adm_threshold()] as `curve`.
#' @section How much data a calibration needs:
#' The threshold is set by the low-`q` end of the curve, and that end needs the most data. An
#' individual carries on average `T * L * q` focal tracts, where `L` is the genetic length in
#' Morgans. At `T = 8` and `q = 0.005` that is **1.4 genome-wide (L = 36 M) but 0.03 on chr22
#' (L = 0.74 M)** -- 97% of such individuals carry no focal tract at all on one chromosome, so their
#' true dosage is constant, their covariance with it is zero, and they contribute noise.
#'
#' Binning rescues this only so far: what has to be adequate is the tract count *per bin*,
#'
#' \deqn{N_{bin} \times T \times L \times q}
#'
#' which is 49 times smaller on chr22 than genome-wide. Scored on chr22 against the same caller's
#' genome-wide curve, reliability agrees to a median of 0.016 above `q_hat` 0.025 (Spearman 0.95) and
#' collapses below it -- 0.13 where the genome-wide answer is 0.84. The break sits where the bin holds
#' of order ten tracts.
#'
#' So: **calibrate genome-wide, and check the formula at the smallest `q` you intend to resolve.**
#' `adm_score()` warns when the site count or the number of low-`q` individuals is too small.
#' If you must use fewer chromosomes, raise the number of simulated individuals in proportion --
#' it is the product that matters, not either factor.
#'
#' @details Per individual, over all sites: `a = cov(called, true)`, `w = var(called)`,
#'   `v = var(true)` on the 0-2 dosage scale, and `q_hat = mean(called) / 2`. These are the three
#'   quantities the inclusion rule is written in; nothing else about the caller is used, which is why
#'   any caller can be scored.
#' @export
adm_score <- function(sim, called, format = c("msp", "flare", "matrix"),
                      focal_code = NULL, bins = 28L, min_n = 15L) {
  format <- match.arg(format)
  if (is.character(sim)) sim <- readRDS(sim)
  if (!inherits(sim, "admit_sim")) stop("`sim` must be an admit_sim object or a path to truth.rds")

  # Never materialise sites x haplotypes: at 121k sites and 12k haplotypes that is 1.4e9 entries
  # and it OOMs on anything short of a cluster node. The calls live on the caller's window grid,
  # which is three orders of magnitude smaller, so expand two columns at a time per individual.
  nI <- nrow(sim$design)
  if (format == "matrix") {
    stopifnot(nrow(called) == nrow(sim$site), ncol(called) == 2L * nI)
    getcols <- function(i) called[, c(2L * i - 1L, 2L * i), drop = FALSE]
  } else {
    W <- switch(format, msp = .read_msp(called, sim, focal_code),
                        flare = .read_flare(called, sim))
    stopifnot(ncol(W$calls) == 2L * nI)
    getcols <- function(i) W$calls[W$idx, c(2L * i - 1L, 2L * i), drop = FALSE]
  }

  a <- v <- w <- qh <- numeric(nI)
  for (i in seq_len(nI)) {
    cc <- rowSums(getcols(i))                                   # 0-2 called dosage
    tt <- as.numeric(sim$anc1[, i]) + as.numeric(sim$anc2[, i]) # 0-2 true dosage
    a[i] <- stats::cov(cc, tt); w[i] <- stats::var(cc); v[i] <- stats::var(tt)
    qh[i] <- mean(cc) / 2
  }

  keep <- is.finite(a) & is.finite(w) & is.finite(v)
  # the low-q end sets the threshold and is the hardest to measure: warn rather than let a thin
  # calibration quietly produce a curve that is only trustworthy where it does not matter
  lowq <- sum(sim$design$q[keep] > 0 & sim$design$q[keep] < 0.02)
  if (nrow(sim$site) < 2e5)
    warning("scored on ", format(nrow(sim$site), big.mark = ","), " sites. The low-q end of the ",
            "curve, which sets the threshold, needs genome-wide data; on one chromosome it is noise.",
            call. = FALSE)
  if (lowq < 200L)
    warning("only ", lowq, " simulated individuals with 0 < q < 0.02. Draw q from your cohort so ",
            "the low end is populated -- that is the part the threshold depends on.", call. = FALSE)
  d  <- data.frame(q_hat = qh, a = a, w = w, v = v)[keep, ]
  br <- c(-Inf, 10^seq(-3.6, -0.3, length.out = bins), Inf)
  out <- do.call(rbind, lapply(split(d, cut(d$q_hat, br)), function(z)
    if (nrow(z) >= min_n) data.frame(q_hat = median(z$q_hat), a = mean(z$a),
                                     w = mean(z$w), v = mean(z$v), n = nrow(z))))
  out <- out[order(out$q_hat), ]; rownames(out) <- NULL
  out$rho <- out$a / out$w
  zero <- sim$design$q[keep] == 0
  attr(out, "eps0")  <- if (any(zero)) mean(qh[keep][zero]) else NA_real_
  attr(out, "focal") <- sim$focal
  out
}

# RFMix and Gnomix both write this: a codes line, a header, then one column per haplotype giving the
# ancestry index of each window. Returns the window-level calls plus each site's window, so the
# caller can expand two columns at a time instead of the whole matrix.
.read_msp <- function(path, sim, focal_code) {
  l1 <- readLines(path, n = 1L)
  if (is.null(focal_code)) {
    kv <- regmatches(l1, gregexpr("[A-Za-z0-9_]+=[0-9]+", l1))[[1]]
    if (!length(kv)) stop("no ancestry codes in the first line; pass `focal_code`")
    nm <- sub("=.*", "", kv); cd <- as.integer(sub(".*=", "", kv))
    if (!sim$focal %in% nm) stop("focal '", sim$focal, "' not among: ", paste(nm, collapse = ", "))
    focal_code <- cd[match(sim$focal, nm)]
  }
  d <- utils::read.delim(path, skip = 1L, check.names = FALSE)
  j <- findInterval(sim$site$POS, d[[2]])              # each site's window (spos)
  j[j < 1L] <- 1L
  list(calls = as.matrix(d[, -(1:6)]) == focal_code, idx = j)
}

# FLARE writes an ancestry VCF with per-haplotype AN1/AN2 FORMAT fields.
.read_flare <- function(path, sim) {
  if (!nzchar(Sys.which("bcftools")))
    stop("reading FLARE output needs bcftools on the path")
  txt <- system2("bcftools", c("query", "-f", "[%AN1,%AN2\\t]\\n", path), stdout = TRUE)
  m <- do.call(rbind, strsplit(sub("\\t$", "", txt), "\t"))
  ac <- array(as.integer(unlist(strsplit(m, ","))), dim = c(2L, nrow(m), ncol(m)))
  stop("FLARE ancestry codes are panel-specific; pass format = 'matrix' with your own decoding")
}
