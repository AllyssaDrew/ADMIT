#' Estimate the background false-call rate from your own cohort
#'
#' The one quantity the inclusion rule needs. `eps0` is what the local-ancestry caller reports for
#' the focal ancestry in people who carry none of it: a truly unadmixed individual still receives
#' invented focal tracts, and that floor does not average away.
#'
#' @section Why the reference group is defined on OTHER ancestries:
#' The obvious recipe -- take people whose *focal* estimate is near zero -- is circular. Those are
#' exactly the individuals the caller got right, so the floor measured on them is biased downwards,
#' and most in a cohort where the focal ancestry is hard to call. The reference group here is instead
#' defined by being near-pure in some *other* ancestry, which does not reference the call being
#' audited. Measured both ways on one cohort the difference was large: the honest recipe kept 38,899
#' individuals and gave 0.449% invented SAS, while the circular one kept 1,956 and would have
#' understated it.
#'
#' @param global_Q Matrix or data frame of global ancestry proportions, one row per individual, one
#'   column per ancestry, columns named. Rows need not sum to exactly 1.
#' @param local_dosage Numeric vector, one value per individual: the genome-wide mean *focal* local
#'   ancestry dosage on a 0-2 scale (or 0-1; see `scale_2`). Same order as `global_Q` rows.
#' @param focal Name of the focal ancestry column in `global_Q`.
#' @param pure_cut An individual joins the reference group when some non-focal ancestry exceeds this.
#' @param scale_2 `TRUE` (default) if `local_dosage` is on a 0-2 diploid scale, so it is halved to a
#'   proportion. Set `FALSE` if you already pass proportions.
#' @param min_ref Error if the reference group is smaller than this.
#'
#' @return An object of class `admit_eps0`: `eps0`, the number and composition of the reference
#'   group, and `q_hat_ref`, their per-individual estimates -- the empirical noise floor, which
#'   [adm_threshold()] uses to report how many unadmixed individuals a cut admits.
#' @examples
#' set.seed(1)
#' n <- 2000
#' Q <- data.frame(AFR = c(runif(200, 0, .8), rep(0, n - 200)))
#' Q$EUR <- 1 - Q$AFR
#' # truly unadmixed individuals still receive invented tracts
#' dose <- 2 * (Q$AFR + abs(rnorm(n, 0, 7e-4)))
#' adm_eps0(Q, dose, focal = "AFR")
#' @export
adm_eps0 <- function(global_Q, local_dosage, focal, pure_cut = 0.99,
                     scale_2 = TRUE, min_ref = 50L) {
  global_Q <- as.data.frame(global_Q)
  if (!focal %in% names(global_Q))
    stop("`focal` must name a column of `global_Q`; have: ", paste(names(global_Q), collapse = ", "))
  if (length(local_dosage) != nrow(global_Q))
    stop("`local_dosage` has ", length(local_dosage), " values for ", nrow(global_Q), " individuals")
  if (pure_cut <= 0.5 || pure_cut >= 1) stop("`pure_cut` must lie in (0.5, 1)")

  other <- setdiff(names(global_Q), focal)
  if (!length(other)) stop("`global_Q` needs at least one ancestry besides `focal`")
  M <- as.matrix(global_Q[, other, drop = FALSE])
  # near-pure in SOME other ancestry -- never defined on the focal estimate itself
  is_ref <- apply(M, 1L, max, na.rm = TRUE) > pure_cut
  is_ref[is.na(is_ref)] <- FALSE
  if (sum(is_ref) < min_ref)
    stop("only ", sum(is_ref), " individuals exceed pure_cut = ", pure_cut,
         " in a non-focal ancestry; lower `pure_cut` or pass `min_ref`")

  qh <- if (scale_2) local_dosage / 2 else local_dosage
  qh_ref <- qh[is_ref]
  qh_ref <- qh_ref[is.finite(qh_ref)]
  which_pure <- other[apply(M[is_ref, , drop = FALSE], 1L, which.max)]

  structure(list(
    eps0       = mean(qh_ref),
    eps0_sd    = sd(qh_ref),
    n_ref      = length(qh_ref),
    n_total    = nrow(global_Q),
    pure_cut   = pure_cut,
    focal      = focal,
    ref_by     = sort(table(which_pure), decreasing = TRUE),
    q_hat_ref  = qh_ref,
    q_hat      = qh
  ), class = "admit_eps0")
}

#' @export
print.admit_eps0 <- function(x, ...) {
  cat(sprintf("Background false-call rate for focal %s\n", x$focal))
  cat(sprintf("  eps0 = %.3g   (sd across the reference group %.3g)\n", x$eps0, x$eps0_sd))
  cat(sprintf("  reference: %s of %s individuals, >%.0f%% in a single non-focal ancestry\n",
              format(x$n_ref, big.mark = ","), format(x$n_total, big.mark = ","), 100 * x$pure_cut))
  cat("  made up of: ",
      paste(sprintf("%s %s", names(x$ref_by), format(as.integer(x$ref_by), big.mark = ",")),
            collapse = ", "), "\n", sep = "")
  cat(sprintf("  implied q* ~ eps0/2 = %.3g  -- a property of the caller, NOT a cut to apply;\n",
              x$eps0 / 2))
  cat(sprintf("  it sits below this estimator's own floor (these individuals average %.3g).\n",
              x$eps0))
  invisible(x)
}
