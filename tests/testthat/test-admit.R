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
  expect_gt(th$detection$cut, 1e-3)              # not the fivefold-too-tight closed-form answer
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
