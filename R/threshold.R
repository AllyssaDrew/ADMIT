#' Reliability of a local-ancestry call at a given global ancestry
#'
#' `r(q) = 2q(1-q) / (2q(1-q) + eps0)`. The numerator is the variance an individual's true local
#' ancestry has across the genome -- the only thing a scan can regress on -- and it vanishes as the
#' individual approaches unadmixed. The denominator adds the caller's noise, which does not vanish.
#' That asymmetry is why an inclusion threshold has to exist.
#'
#' @param q Global focal ancestry proportion(s).
#' @param eps0 Background false-call rate, from [adm_eps0()].
#' @param lambda Calling gain, `cov(called, true)/var(true)`. Hard-call misclassification makes this
#'   below 1 (0.934 measured for RFMix). A gain that does not vary with `q` costs no detection power
#'   -- the non-centrality parameter is invariant to rescaling the called dosage -- so it is carried
#'   here only so that effect-size attenuation is right.
#' @return Reliability on \[0, 1\].
#' @export
adm_reliability <- function(q, eps0, lambda = 1) {
  v <- 2 * q * (1 - q)
  (lambda * v) / (lambda^2 * v + eps0)
}

# Per-individual signal and variance contributions under the closed-form error model.
#
# `q` arriving here is the ESTIMATE, and the estimate is biased upward by the floor: a truly
# unadmixed individual still reads eps0. Feeding the estimate straight into 2q(1-q) therefore counts
# the caller's noise as signal, and the consequence is not subtle -- on a cohort of 52,639 where
# 39,069 individuals are near-pure in another ancestry, it credits each of them with reliability
# 0.66 when their true reliability is 0, and the optimum lands an order of magnitude too low.
# So the floor is removed first: E[q_hat | q] ~ lambda*q + eps0 for small q, inverted and truncated.
.aw <- function(q_hat, eps0, lambda) {
  q <- pmax((q_hat - eps0) / lambda, 0)
  v <- 2 * q * (1 - q)
  list(a = lambda * v, w = lambda^2 * v + eps0, v = v, q = q)
}

#' Efficiency of an inclusion window
#'
#' The non-centrality parameter of the scan is `(sum a)^2 / sum w` over the included individuals,
#' where `a` is the signal each contributes and `w` the variance. Efficiency divides that by the
#' value an error-free caller would reach on the whole cohort, so 1 is the ceiling.
#'
#' @param q_hat Estimated global focal ancestry, one value per individual -- the ESTIMATE, not the
#'   truth. Selecting on the estimate is not a concession: it is what you can do, and it is also
#'   better, because the quantity worth selecting on is reliability and the estimate is partly a
#'   readout of it while the truth is not.
#' @param eps0 Background false-call rate, from [adm_eps0()].
#' @param cuts Candidate lower cuts. Default is a log grid spanning the cohort.
#' @param lambda Calling gain; see [adm_reliability()].
#' @param upper Upper cut, applied symmetrically as `q_hat < upper`.
#' @param curve A measured reliability curve from [adm_curve()] or [adm_rescale()]. Supplying one
#'   replaces the closed form, which gets the scaling right and the level wrong.
#' @return Data frame: `cut`, `n` retained, `eff`, and `R`, the attenuation factor `sum a / sum w`.
#' @export
adm_efficiency <- function(q_hat, eps0 = NULL, cuts = NULL, lambda = 1, upper = 1, curve = NULL) {
  if (is.null(curve) && is.null(eps0)) stop('give either `eps0` (closed form) or `curve` (measured)')
  q_hat <- q_hat[is.finite(q_hat)]
  q <- pmin(q_hat, 1 - q_hat)                       # fold: 0.98 focal is as informative as 0.02
  if (is.null(cuts)) cuts <- 10^seq(log10(max(min(q[q > 0]), 1e-5)), log10(0.2), length.out = 120)
  x <- if (is.null(curve)) .aw(q, eps0, lambda) else .aw_curve(q, curve)
  Imax <- sum(x$v)                                  # error-free, whole cohort
  out <- lapply(cuts, function(c0) {
    s <- q > c0 & q_hat < upper                     # ... but the CUT is applied to the estimate
    if (!any(s)) return(data.frame(cut = c0, n = 0L, eff = 0, R = NA_real_))
    A <- sum(x$a[s]); W <- sum(x$w[s])
    data.frame(cut = c0, n = sum(s), eff = A^2 / W / Imax, R = A / W)
  })
  do.call(rbind, out)
}

#' Choose the inclusion window
#'
#' Returns the cut that maximises detection power and the cut that minimises the mean squared error
#' of the effect estimate. They differ because loosening the window moves them oppositely: the
#' variance term improves as more individuals enter, while attenuation worsens. The interval between
#' them is the defensible range.
#'
#' @param q_hat Estimated global focal ancestry per individual.
#' @param eps0 Either a number or an `admit_eps0` object from [adm_eps0()]. Passing the object also
#'   lets the result report how many unadmixed individuals each cut admits.
#' @param lambda Calling gain; see [adm_reliability()].
#' @param ncp_target Fixes the unknown ratio of residual to effect variance by declaring what
#'   non-centrality the widest window would carry. Only the estimation optimum depends on it; the
#'   detection optimum does not.
#' @param cuts Candidate cuts; default a log grid.
#' @param curve A measured reliability curve from [adm_curve()] or [adm_rescale()], or the string
#'   `"AFR"` / `"SAS"`. Strongly preferred over the closed form; see [adm_curve()].
#' @param upper Upper cut.
#' @return An object of class `admit_threshold`.
#' @examples
#' set.seed(1)
#' q <- c(runif(300, 0, .9), abs(rnorm(5000, 0, 1.3e-3)))
#' adm_threshold(q, eps0 = 1.4e-3)
#' @export
adm_threshold <- function(q_hat, eps0 = NULL, curve = NULL, lambda = 1, ncp_target = 25,
                          cuts = NULL, upper = 1) {
  e <- eps0
  floor_q <- NULL
  if (inherits(eps0, "admit_eps0")) { floor_q <- eps0$q_hat_ref; e <- eps0$eps0 }
  if (is.character(curve)) curve <- adm_curve(curve)
  if (is.null(curve)) {
    if (!is.numeric(e) || length(e) != 1L || !is.finite(e) || e <= 0)
      stop("give `curve`, or an `eps0` that is a positive number or an `admit_eps0` object")
    warning("the one-parameter closed form reproduces how the optimum SCALES with eps0 but not ",
            "its level: against known truth it sits about 5x too tight and overstates efficiency ",
            "by ~0.18. Prefer curve = 'AFR' / 'SAS', or calibrate your own.", call. = FALSE)
  }
  tab <- adm_efficiency(q_hat, e, cuts = cuts, lambda = lambda, upper = upper, curve = curve)
  tab <- tab[tab$n > 0, , drop = FALSE]
  if (!nrow(tab)) stop("no cut retains anyone")

  # detection: maximise the NCP. estimation: minimise MSE/gamma^2 = (1-R)^2 + (sigma^2/gamma^2)/W,
  # with sigma^2/gamma^2 pinned by `ncp_target` at the widest window.
  i_det <- which.max(tab$eff)

  qf <- pmin(q_hat, 1 - q_hat); qf <- qf[is.finite(qf)]
  x  <- if (is.null(curve)) .aw(qf, e, lambda) else .aw_curve(qf, curve)
  Av <- vapply(tab$cut, function(c0) sum(x$a[qf > c0]), 0)
  Wv <- vapply(tab$cut, function(c0) sum(x$w[qf > c0]), 0)
  # NCP = (gamma^2/sigma^2) * A^2/W. Pin the unknown ratio by declaring the NCP the widest window
  # carries, then sigma^2/gamma^2 = A_wide^2 / (ncp_target * W_wide). Only the estimation optimum
  # depends on this; the detection optimum does not.
  j <- which.max(tab$n)
  s2g2 <- Av[j]^2 / (ncp_target * Wv[j])
  tab$mse <- (1 - tab$R)^2 + s2g2 / Wv
  i_est <- which.min(tab$mse)

  # how many truly-unadmixed individuals each cut admits -- the number that stops q* being used
  if (!is.null(floor_q) && length(floor_q)) {
    tab$admits_unadmixed <- vapply(tab$cut, function(c0) mean(floor_q > c0), 0)
  }

  structure(list(
    table = tab, eps0 = e, lambda = lambda, ncp_target = ncp_target,
    detection = tab[i_det, ], estimation = tab[i_est, ],
    qstar = if (is.null(e)) NA_real_ else e / 2, n_total = sum(is.finite(q_hat)),
    mode = if (is.null(curve)) "one-parameter closed form" else "measured curve",
    has_floor = !is.null(floor_q), floor_q = floor_q
  ), class = "admit_threshold")
}

#' @export
print.admit_threshold <- function(x, ...) {
  d <- x$detection; s <- x$estimation
  lo <- min(d$cut, s$cut); hi <- max(d$cut, s$cut)
  cat(sprintf("Inclusion window  [%s]\n\n", x$mode))
  cat(sprintf("  detection optimum    %-9.3g  efficiency %.3f   keeps %s\n",
              d$cut, d$eff, format(d$n, big.mark = ",")))
  cat(sprintf("  estimation optimum   %-9.3g  efficiency %.3f   keeps %s\n",
              s$cut, s$eff, format(s$n, big.mark = ",")))
  cat(sprintf("\n  recommended window   %.3g to %.3g   (anything inside is within reach of both)\n",
              lo, hi))
  if (x$has_floor) {
    f <- function(c0) mean(x$floor_q > c0)
    cat(sprintf("\n  unadmixed individuals admitted: %.1f%% at the detection cut\n", 100 * f(d$cut)))
    if (is.finite(x$qstar))
      cat(sprintf("  q* = %.3g admits %.1f%% of them -- it tracks the caller's error, it is not a cut.\n",
                  x$qstar, 100 * f(x$qstar)))
  }
  invisible(x)
}
