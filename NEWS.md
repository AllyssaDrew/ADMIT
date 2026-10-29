# ADMIT 0.1.1

* `adm_rescale()` on the shipped South Asian curve scaled from the African curve's `eps0`: `adm_curve()` dropped
  the curve's ancestry label, so the default fell back to AFR (1.4e-3) instead of the South Asian simulation's
  2.46e-3, misplacing the rescaled SAS curve by a factor of 1.76. Fixed, and the stored reference values are now
  those of the simulations the curves were measured on: AFR 1.33e-3 (was 1.4e-3), SAS 2.46e-3 (was 4.4e-3, the
  real UK Biobank value). A curve that is not a shipped one now requires `eps0_from`.
* Documentation: the shipped curves come from 3,600 genome-wide mosaics per ancestry with the low end
  oversampled, not 16,000 mosaics drawn from the UK Biobank distribution; `adm_curve()` no longer calls
  rescaling untested.

# ADMIT 0.1.0

First release.

* `adm_eps0()` estimates the caller's background false-call rate from the user's own cohort, with no
  simulation. The reference group is defined by being near-pure in a *non-focal* ancestry, because
  defining it on the focal estimate selects the individuals the caller got right.
* `adm_threshold()` returns the detection- and estimation-optimal inclusion cuts, the efficiency
  cost of any other cut, and how many unadmixed individuals each cut admits.
* `adm_curve()` ships reliability curves measured against held-out-founder truth for RFMix 2 on the
  5-way HGDP+1kGP reference. `adm_rescale()` moves a curve to another caller by the ratio of their
  `eps0`; validated against FLARE, where it costs 0.001 of efficiency and skipping it costs 0.016.
* `adm_holdout()`, `adm_simulate()` and `adm_score()` calibrate a curve for any caller. The package
  generates the query and scores the output; it does not run the caller.
* The one-parameter closed form is available but warns: it reproduces how the optimum scales with
  `eps0` and not its level, sitting about twofold too tight against known truth.
