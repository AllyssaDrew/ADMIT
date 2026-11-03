test_that("adm_eps0 recovers a planted floor and refuses a circular reference group", {
  set.seed(7); n <- 4000
  Q <- data.frame(AFR = c(runif(300, 0, .8), rep(0, n - 300)))
  Q$EUR <- 1 - Q$AFR
  truth <- 8e-4
  dose <- 2 * (Q$AFR + abs(rnorm(n, 0, truth * sqrt(pi/2))))
  e <- adm_eps0(Q, dose, focal = "AFR")
  expect_s3_class(e, "admit_eps0")
  expect_equal(e$eps0, truth, tolerance = 0.25)
  expect_gt(e$n_ref, 3000)                       # the near-pure EUR individuals
  expect_error(adm_eps0(Q[, "AFR", drop = FALSE], dose, focal = "AFR"), "besides")
})

test_that("reliability is 0 for an unadmixed individual and approaches 1 when admixed", {
  expect_equal(adm_reliability(0, eps0 = 1e-3), 0)
  expect_gt(adm_reliability(0.3, eps0 = 1e-3), 0.99)
  expect_lt(adm_reliability(1e-4, eps0 = 1e-3), 0.2)
  # Reliability depends on q essentially through q/eps0 -- the property adm_rescale() relies on.
  # The dependence is NOT exact: r = 2q(1-q)/(2q(1-q)+eps0), and the (1-q) factor breaks scale
  # invariance at O(q). It is exact in the limit and good to ~1e-3 over the range that matters,
  # which is why the empirical collapse of the RFMix and FLARE curves is 0.005 and not 0.
  expect_equal(adm_reliability(2e-4, 2e-4), adm_reliability(4e-4, 4e-4), tolerance = 1e-4)
  expect_equal(adm_reliability(2e-3, 2e-3), adm_reliability(4e-3, 4e-3), tolerance = 1e-3)
  expect_gt(abs(adm_reliability(0.05, 0.05) - adm_reliability(0.10, 0.10)), 1e-3)  # degrades at large q
})

test_that("the shipped curve reproduces the published UK Biobank window", {
  set.seed(1)
  th <- adm_threshold(c(runif(3000, 0, .9), abs(rnorm(49000, 0, 1.1e-3))), curve = "AFR")
  expect_s3_class(th, "admit_threshold")
  expect_equal(th$mode, "measured curve")
  expect_gt(th$detection$cut, 1e-3)              # the measured curve, not the closed form (about twofold too tight)
  expect_lt(th$detection$cut, 2e-2)
  expect_true(th$estimation$cut >= th$detection$cut)   # estimation always the looser of the two
})

test_that("the closed form warns, because it gets the level wrong", {
  expect_warning(adm_threshold(c(runif(500, 0, .5), abs(rnorm(5000, 0, 1e-3))), eps0 = 1.4e-3),
                 "SCALES with eps0 but not")
})

test_that("adm_rescale moves the axis by exactly the eps0 ratio", {
  cur <- adm_curve("AFR")
  out <- adm_rescale(cur, eps0_from = 1.33e-3, eps0_to = 3.19e-3)
  expect_equal(out$q_hat, cur$q_hat * (3.19e-3 / 1.33e-3))
  expect_equal(out$a, cur$a)                     # only the axis moves
  expect_error(adm_rescale(cur, 1.33e-3, -1), "positive")
})

test_that("efficiency is bounded and monotone in the right places", {
  set.seed(3); q <- c(runif(500, 0, .5), abs(rnorm(5000, 0, 1e-3)))
  tab <- adm_efficiency(q, curve = adm_curve("AFR"))
  expect_true(all(tab$eff >= 0 & tab$eff <= 1))
  expect_true(all(diff(tab$n) <= 0))             # a stricter cut never keeps more people
})

test_that("adm_rescale scales each shipped curve from its own simulation eps0", {
  k <- attr(ADMIT:::adm_curves, "eps0")
  k <- k[!is.na(k)]                      # the American curves carry none on purpose; see their own test
  for (key in names(k)) {
    p <- strsplit(key, "|", fixed = TRUE)[[1]]                       # caller | focal | background
    cur <- adm_curve(p[2], caller = p[1], background = p[3])
    expect_equal(attr(adm_rescale(cur, eps0_to = k[[key]]), "rescaled_by"), 1)
  }
  expect_equal(unname(k[c("RFMix|AFR|European", "RFMix|SAS|European")]), c(1.33e-3, 2.46e-3), tolerance = 0.01)
})

test_that("a curve belongs to a pair of ancestries, and rescaling warns across a wide gap", {
  eu <- adm_curve("SAS", caller = "FLARE", background = "European")
  af <- adm_curve("SAS", caller = "FLARE", background = "African")
  r <- function(cur, q) approx(log(cur$q_hat), cur$rho, xout = log(q), rule = 2)$y
  expect_gt(r(af, 0.01), 0.9)        # South Asian on an African background is recovered
  expect_lt(r(eu, 0.01), 0.2)        # the same ancestry on a European background is not
  expect_warning(adm_rescale(adm_curve("AFR", caller = "FLARE"), eps0_to = 6e-3), "rescaling by")
  expect_silent(adm_rescale(adm_curve("AFR", caller = "FLARE"), eps0_to = 1e-3))
})

test_that("adm_eps0(background=) restricts the reference group to that ancestry", {
  set.seed(3); n <- 3000
  Q <- data.frame(AFR = c(runif(300, 0, .8), rep(0, n - 300)))
  Q$EUR <- c(1 - Q$AFR[1:300], rep(1, 1350), rep(0, n - 1650))
  Q$SAS <- 1 - Q$AFR - Q$EUR
  dose <- 2 * (Q$AFR + ifelse(Q$SAS > 0.99, 6e-3, 6e-4) * abs(rnorm(n)))
  any <- adm_eps0(Q, dose, focal = "AFR")
  eur <- adm_eps0(Q, dose, focal = "AFR", background = "EUR")
  sas <- adm_eps0(Q, dose, focal = "AFR", background = "SAS")
  expect_equal(eur$n_ref + sas$n_ref, any$n_ref)       # the two backgrounds partition the reference group
  expect_gt(sas$eps0, 5 * eur$eps0)                     # and they measure different noise floors
  expect_identical(eur$background, "EUR")
  expect_error(adm_eps0(Q, dose, focal = "AFR", background = "AFR"), "non-focal")
})

test_that("adm_threshold reports a plateau that contains the optimum and ends at a real cut", {
  set.seed(11)
  q <- c(runif(400, 0, .9), abs(rnorm(6000, 0, 2e-3)))
  th <- adm_threshold(q, curve = adm_curve("AFR"))
  expect_length(th$plateau, 2L)
  expect_lte(th$plateau[1], th$detection$cut)
  expect_gte(th$plateau[2], th$detection$cut)
  expect_equal(th$conservative$cut, th$plateau[2])
  expect_gte(th$conservative$eff, th$detection$eff - th$plateau_tol)   # free, by definition of the plateau
  expect_lte(th$conservative$n, th$detection$n)                        # and it is the tighter end
  wide <- adm_threshold(q, curve = adm_curve("AFR"), plateau_tol = 0.2)
  expect_gte(diff(wide$plateau), diff(th$plateau))                     # a looser tolerance cannot narrow it
})

test_that("the American curves are shipped as measured and cannot be rescaled", {
  eu <- adm_curve("AMR", caller = "FLARE", background = "European")
  af <- adm_curve("AMR", caller = "FLARE", background = "African")
  expect_gt(nrow(eu), 20); expect_gt(nrow(af), 20)
  r <- function(cur, q) approx(log(cur$q_hat), cur$rho, xout = log(q), rule = 2)$y
  expect_lt(r(eu, 0.002), 0.1)      # the window that was in use, 0.00076, sits where reliability is ~0
  expect_gt(r(eu, 0.1), 0.8)
  expect_gt(r(af, 0.03), r(eu, 0.03))        # American is easier to call against an African background
  expect_error(adm_rescale(eu, eps0_to = 1.33e-4), "used as measured")
  expect_error(adm_rescale(af, eps0_to = 1.33e-4), "used as measured")
})

test_that("a plateau that reaches the grid's edge is marked censored", {
  set.seed(4); q <- c(runif(400, 0, .9), abs(rnorm(6000, 0, 2e-3)))
  wide <- adm_threshold(q, curve = adm_curve("AFR"))
  expect_false(any(wide$plateau_censored))                       # the default grid is wide enough here
  narrow <- adm_threshold(q, curve = adm_curve("AFR"), cuts = 10^seq(-2.6, -2.2, length.out = 20))
  expect_true(all(narrow$plateau_censored))                      # this one cannot see either boundary
  expect_output(print(narrow), "search grid")
})
