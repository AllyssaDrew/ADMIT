#' Reliability curves measured against known truth
#'
#' The closed form in [adm_reliability()] has one free parameter and gets the *scaling* right --
#' the optimal cut moves in proportion to `eps0` -- but not the *level*. Measured against a
#' held-out-founder simulation on the same cohort it is too optimistic at the low end (reliability
#' 0.58 where the truth is 0.14) and places the optimum about fivefold too tight. So the default is
#' not a formula, it is a measured curve.
#'
#' @param focal `"AFR"` or `"SAS"`.
#' @return Data frame of `q_hat`, `a`, `w`, `v`, `rho`, `n`.
#' @details Shipped from RFMix 2 against the 5-way HGDP+1kGP reference (AFR/EAS/EUR/NAM/SAS),
#'   16,000 held-out-founder mosaics, T = 8 generations, q drawn from the real UK Biobank
#'   distribution. **A different caller or reference panel is a different curve.** Whether a
#'   shipped curve can simply be rescaled by the ratio of the two `eps0` values is untested; until
#'   it is, treat the shipped curve as a starting point and calibrate your own if your setup
#'   differs.
#' @export
adm_curve <- function(focal = c("AFR", "SAS")) {
  focal <- match.arg(focal)
  d <- adm_curves[adm_curves$focal == focal, , drop = FALSE]
  rownames(d) <- NULL
  d
}

# interpolate a measured curve onto the individuals' q_hat; flat extrapolation at both ends
.aw_curve <- function(q_hat, curve) {
  f <- function(y) approx(curve$q_hat, y, xout = q_hat, rule = 2)$y
  list(a = f(curve$a), w = f(curve$w), v = f(curve$v), q = q_hat)
}

#' Rescale a measured curve to a different caller
#'
#' Reliability depends on the global ancestry only through `q / eps0`, so a curve measured for one
#' caller can in principle be moved to another by scaling its `q_hat` axis by the ratio of their
#' background false-call rates. Whether that holds for *measured* curves rather than the closed form
#' is an empirical question, and it was tested.
#'
#' @section What the test showed:
#' RFMix 2 and FLARE were run on the same 16,000 held-out-founder mosaics, where they differ 2.4-fold
#' in `eps0` (1.33e-3 against 3.19e-3). Plotted against the scaled axis `q/eps0` the two measured
#' reliability curves agree to a median of **0.005** across 0.38 < q/eps0 < 319. Operationally, on a
#' real 52,639-individual cohort: a FLARE user picking their window from FLARE's own curve attains
#' efficiency 0.859; picking it from the shipped RFMix curve rescaled attains **0.858**, and from the
#' shipped curve used unrescaled, 0.843. So rescaling costs 0.001 and skipping it costs 0.016.
#'
#' Two things that follow. Rescaling is a supported route, so a caller we have not measured does not
#' force a full calibration. And the choice is forgiving: the efficiency curve is a plateau near its
#' optimum, so being twofold off on the cut costs under 0.02. Report the window to one significant
#' figure, not two.
#'
#' @param curve A curve from [adm_curve()], or any data frame with `q_hat`, `a`, `w`, `v`.
#' @param eps0_from `eps0` of the caller the curve was measured on. Defaults to the shipped value
#'   when `curve` came from [adm_curve()].
#' @param eps0_to `eps0` of your caller, from [adm_eps0()] -- a number or an `admit_eps0` object.
#' @return The curve with its `q_hat` axis scaled by `eps0_to / eps0_from`.
#' @examples
#' cur <- adm_curve("AFR")
#' head(adm_rescale(cur, eps0_from = 1.33e-3, eps0_to = 3.19e-3), 3)
#' @export
adm_rescale <- function(curve, eps0_from = NULL, eps0_to) {
  if (is.character(curve)) curve <- adm_curve(curve)
  if (inherits(eps0_to, "admit_eps0")) eps0_to <- eps0_to$eps0
  if (is.null(eps0_from)) eps0_from <- attr(adm_curves, "eps0")[[attr(curve, "focal") %||% "AFR"]]
  if (!is.numeric(eps0_from) || !is.numeric(eps0_to) || eps0_from <= 0 || eps0_to <= 0)
    stop("`eps0_from` and `eps0_to` must be positive numbers")
  curve$q_hat <- curve$q_hat * (eps0_to / eps0_from)
  attr(curve, "rescaled_by") <- eps0_to / eps0_from
  curve
}

`%||%` <- function(a, b) if (is.null(a)) b else a
