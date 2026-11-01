#' Reliability curves measured against known truth
#'
#' The closed form in [adm_reliability()] has one free parameter and gets the *scaling* right --
#' the optimal cut moves in proportion to `eps0` -- but not the *level*. Measured against a
#' held-out-founder simulation on the same cohort it is too optimistic at the low end (reliability
#' 0.59 where the measured value is 0.15) and places the optimum about twofold too tight. So the default is
#' not a formula, it is a measured curve.
#'
#' @param focal `"AFR"` or `"SAS"`.
#' @param caller `"RFMix"` (default) or `"FLARE"`.
#' @param background the ancestry the rest of the genome comes from: `"European"` (default) or,
#'   for South Asian ancestry called by FLARE, `"African"`.
#' @return Data frame of `q_hat`, `a`, `w`, `v`, `rho`, `n`.
#' @details Two callers are shipped.
#'   `RFMix` is RFMix 2 against the 5-way 446-sample HGDP+1kGP reference used in UK Biobank:
#'   genome-wide held-out-founder mosaics, 3,600 per focal ancestry, single pulse T = 8 generations.
#'   `FLARE` is FLARE 0.6.0 against a 3,201-sample 5-way reference matching the All of Us setup:
#'   16,000 mosaics, founders held out in five rotating folds with their 1000 Genomes relatives removed,
#'   q covering the whole 0-1 range, validated against the real 52,639-participant UK Biobank cohort
#'   run through the same FLARE setup.
#'
#'   **A curve belongs to a pair of ancestries, not to one.** Reliability depends on what the rest of the
#'   genome is, because the caller has to tell the focal ancestry apart from it. South Asian ancestry on an
#'   African background is recovered almost perfectly (reliability 0.98 at q = 0.01); the same ancestry on a
#'   European background is barely recovered at all (0.09), because South Asian and European haplotypes are
#'   far less differentiated. African ancestry on a European background sits at 0.68. In a cohort whose
#'   participants have mixed backgrounds a single curve averages over groups that differ this much, so pick
#'   the background your cohort actually has, or stratify.
#'
#'   Each curve carries the `eps0` of the simulation it was measured on, computed the way [adm_eps0()]
#'   computes it on real data, which is what [adm_rescale()] scales from.
#' @export
adm_curve <- function(focal = c("AFR", "SAS"), caller = c("RFMix", "FLARE"),
                      background = c("European", "African")) {
  focal <- match.arg(focal); caller <- match.arg(caller); background <- match.arg(background)
  key <- paste(caller, focal, background, sep = "|")
  d <- adm_curves[adm_curves$focal == focal & adm_curves$caller == caller &
                    adm_curves$background == background, , drop = FALSE]
  if (!nrow(d)) {
    have <- unique(paste(adm_curves$caller, adm_curves$focal, adm_curves$background, sep = "|"))
    stop("no shipped curve for ", key, "; available: ", paste(have, collapse = ", "))
  }
  rownames(d) <- NULL
  attr(d, "curve_key") <- key                   # adm_rescale() reads the curve's own eps0 from this
  attr(d, "focal") <- focal
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
#' RFMix 2 and FLARE were run on the same genome-wide held-out-founder mosaics (3,600 per focal ancestry), where they differ 2.4-fold
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
#' head(adm_rescale(cur, eps0_from = 1.33e-3, eps0_to = 1.6e-3), 3)
#' @export
adm_rescale <- function(curve, eps0_from = NULL, eps0_to) {
  if (is.character(curve)) curve <- adm_curve(curve)
  if (inherits(eps0_to, "admit_eps0")) eps0_to <- eps0_to$eps0
  if (is.null(eps0_from)) {
    k <- attr(curve, "curve_key")
    if (is.null(k) || !k %in% names(attr(adm_curves, "eps0")))
      stop("`eps0_from` is needed: this curve is not one of the shipped curves")
    eps0_from <- attr(adm_curves, "eps0")[[k]]
  }
  if (!is.numeric(eps0_from) || !is.numeric(eps0_to) || eps0_from <= 0 || eps0_to <= 0)
    stop("`eps0_from` and `eps0_to` must be positive numbers")
  r <- eps0_to / eps0_from
  if (r > 2 || r < 0.5)
    warning("rescaling by ", signif(r, 3), "-fold. A rescale slides the curve along q; it cannot change ",
            "how much of the true dosage the caller recovers, so over a gap this wide it will misstate ",
            "reliability where q is large. Measured example: the AFR curve rescaled to the eps0 of South ",
            "Asian ancestry on a European background reports 0.94 at q = 0.6 where the measured value is ",
            "0.81, because that pairing plateaus at 0.84 and no rescale can lower a plateau. ",
            "Calibrate directly (adm_simulate) if the window matters.", call. = FALSE)
  curve$q_hat <- curve$q_hat * r
  attr(curve, "rescaled_by") <- r
  curve
}

`%||%` <- function(a, b) if (is.null(a)) b else a
