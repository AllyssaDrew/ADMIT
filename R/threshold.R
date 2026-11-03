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
#'   truth. A cut on the estimate attains the same efficiency as a cut on the unobservable truth.
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
#' @param plateau_tol Cuts whose efficiency is within this of the optimum form the plateau: the set the power
#'   criterion cannot tell apart. Default 0.02, the scale over which the efficiency curve is flat in practice.
#' @param upper Upper cut.
#' @return An object of class `admit_threshold`: the full table, the detection and estimation optima, and
#'   `plateau` / `conservative` -- the range of cuts power cannot separate and the tightest of them -- and
#'   `plateau_censored`, naming any end of that range that is the search grid's edge rather than a measured
#'   boundary. On a cohort whose participants reach very low `q_hat`, the lower end usually is.
#' @examples
#' set.seed(1)
#' q <- c(runif(300, 0, .9), abs(rnorm(5000, 0, 1.3e-3)))
#' adm_threshold(q, eps0 = 1.4e-3)
#' @export
adm_threshold <- function(q_hat, eps0 = NULL, curve = NULL, lambda = 1, ncp_target = 25,
                          cuts = NULL, upper = 1, plateau_tol = 0.02) {
  e <- eps0
  floor_q <- NULL
  if (inherits(eps0, "admit_eps0")) { floor_q <- eps0$q_hat_ref; e <- eps0$eps0 }
  if (is.character(curve)) curve <- adm_curve(curve)
  if (is.null(curve)) {
    if (!is.numeric(e) || length(e) != 1L || !is.finite(e) || e <= 0)
      stop("give `curve`, or an `eps0` that is a positive number or an `admit_eps0` object")
    warning("the one-parameter closed form reproduces how the optimum SCALES with eps0 but not ",
            "its level: against known truth it sits about 2x too tight and overstates efficiency ",
            "by ~0.13. Prefer curve = 'AFR' / 'SAS', or calibrate your own.", call. = FALSE)
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

  # The plateau. Reporting one cut invites the reader to think it is special; over a wide range of cuts the
  # efficiency barely moves, and hiding that hides the fact that the choice is still open. The plateau is every
  # cut whose efficiency is within `plateau_tol` of the optimum -- the set power cannot tell apart. Its upper
  # end is the tightest cut power still treats as free, which matters because the individuals it drops are the
  # ones whose calls carry the most error; whether to take it is a decision this function does not make.
  keep <- tab$eff >= tab$eff[i_det] - plateau_tol
  plat <- range(tab$cut[keep])
  # A plateau that reaches an end of the search grid is CENSORED there: the number printed would be the grid's
  # edge, not a boundary anyone measured, and it looks like a measurement. Record which ends are real.
  plat_censored <- c(lo = isTRUE(all.equal(plat[1], min(tab$cut))),
                     hi = isTRUE(all.equal(plat[2], max(tab$cut))))
  i_cons <- which(tab$cut == max(tab$cut[keep]))[1]

  structure(list(
    table = tab, eps0 = e, lambda = lambda, ncp_target = ncp_target,
    detection = tab[i_det, ], estimation = tab[i_est, ],
    plateau = plat, plateau_tol = plateau_tol, plateau_censored = plat_censored,
    conservative = tab[i_cons, ],
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
  if (!is.null(x$plateau)) {
    cn <- x$conservative
    ce <- x$plateau_censored %||% c(lo = FALSE, hi = FALSE)
    # A censored end is NOT printed as a number with a qualifier -- "<=1e-05 to 0.032" still reads as a range
    # someone measured, and gets quoted that way. Say what was actually established instead.
    if (!any(ce)) {
      cat(sprintf("\n  plateau              %.3g to %.3g   (efficiency within %.2f of the optimum)\n",
                  x$plateau[1], x$plateau[2], x$plateau_tol))
    } else {
      cat(sprintf("\n  plateau              efficiency stays within %.2f of the optimum for every cut from\n",
                  x$plateau_tol))
      cat(sprintf("                       %.3g, the smallest evaluated, up to %.3g%s\n", x$plateau[1],
                  x$plateau[2], if (ce[["hi"]]) ", the largest evaluated" else ""))
      cat(sprintf("                       no %s bound within the evaluated range -- pass `cuts` to look further\n",
                  if (all(ce)) "lower or upper" else if (ce[["lo"]]) "lower" else "upper"))
    }
    cat(sprintf("  conservative end     %-9.3g  efficiency %.3f   keeps %s\n",
                cn$cut, cn$eff, format(cn$n, big.mark = ",")))
    cat("  Power cannot separate cuts inside the plateau. The conservative end costs at most\n")
    cat(sprintf("  %.2f of efficiency and drops %s individuals whose calls are least reliable.\n",
                x$plateau_tol, format(d$n - cn$n, big.mark = ",")))
  }
  if (x$has_floor) {
    f <- function(c0) mean(x$floor_q > c0)
    cat(sprintf("\n  unadmixed individuals admitted: %.1f%% at the detection cut\n", 100 * f(d$cut)))
    if (is.finite(x$qstar))
      cat(sprintf("  q* = %.3g admits %.1f%% of them -- it tracks the caller's error, it is not a cut.\n",
                  x$qstar, 100 * f(x$qstar)))
  }
  invisible(x)
}
