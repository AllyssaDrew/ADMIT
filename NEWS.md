# ADMIT 0.4.2

* `adm_efficiency()` accepts a curve by name (`curve = "AFR"` or `"SAS"`), as `adm_threshold()` does and as the
  vignette shows. Until now it needed the data frame from `adm_curve()` and stopped with an error on a name.
* README and vignette now show the results of the curves shipped since 0.4.0. On UK Biobank African ancestry the
  window is 0.0020 to 0.0024 (efficiency 0.893), not 0.0038 to 0.0057 (0.840); for South Asian ancestry 0.0084 to
  0.025 (0.494). The seven shipped curves are listed with their caller, reference panel and number of genomes.
* The section on moving a curve to another caller is rewritten. The two shipped curves for African ancestry were
  measured with different reference panels and no longer coincide once rescaled (median difference 0.06, larger at
  low global ancestry). Rescaling still recovers most of the efficiency: 0.908 against 0.914 from FLARE's own curve
  and 0.855 without rescaling. The text no longer says that the rescaled curves agree to 0.005.
* No change to any computed result.

# ADMIT 0.4.1

* `adm_efficiency()` and `adm_threshold()` take `call_var`, the variance of each individual's called dosage along
  the genome. The curves assume that individuals with the same `q_hat` vary along the genome alike. Above a global
  ancestry of 30% that overstates the variance of real calls (by 21% to 38% for RFMix in UK Biobank), mostly
  because a person with one parent of each ancestry carries one focal copy at every position and so informs no
  scan. With `call_var`, each individual above `call_var_from` (default 0.3) is weighted by the variance observed
  in place of the variance assumed; reliability is unchanged. On UK Biobank African ancestry the optimum does not
  move; the predicted efficiency and attenuation fall. Without `call_var` every result is as in 0.4.0.

# ADMIT 0.4.0

* The two RFMix curves (`adm_curve("AFR")`, `adm_curve("SAS")`) are replaced. They now come from mosaics of gnomAD
  individuals who are absent from the reference panel, on a north-west European background, called with the full
  446-sample reference, and reliability is measured over RFMix windows, the unit an admixture scan tests. The
  previous pair was measured with the founders removed from the reference, which raised the background false-call
  rate to 2.3 times the value seen in UK Biobank, and over genotyped sites. On the UK Biobank distribution of
  African ancestry the detection optimum moves from 0.38% to 0.20% and the efficiency there from 0.840 to 0.893;
  for South Asian ancestry from 0.54% to 0.84% (0.446 to 0.494). Reference `eps0`: AFR 5.12e-4, SAS 3.31e-3.
* A planted-truth simulation (effect on true local ancestry, scan on calls, false calls placed as in real data)
  agrees with `adm_efficiency()` at eight cuts: detection within 0.009, attenuation within 0.008.

# ADMIT 0.3.2

* A plateau that reaches an end of the search grid is now marked as censored there (`plateau_censored`), and
  is no longer printed as a range. The grid's lower bound is `max(min(q_hat > 0), 1e-5)`, so on a cohort whose
  participants reach very low global ancestry the plateau's lower end is usually the grid's edge; printing
  "1e-05 to 0.032" reads as a range someone measured and gets quoted that way. The censored form states what
  was established instead -- efficiency stays within the tolerance for every cut from the smallest evaluated
  up to the upper end, with no lower bound inside the evaluated range.
# ADMIT 0.3.1

* Measured American-ancestry curves (`adm_curve("AMR", caller = "FLARE")`), on a European or African
  background, from an 8,000-mosaic arm. They replace what was being done before, which was to rescale the
  African curve 8.7-fold -- and the measurement shows why that could not work: the calling gain for American
  ancestry is 0.895 against 0.979 for African, so the two differ in more than a scale factor, and the rescaled
  curve put the window where measured reliability is 0.
* **The American curves carry no `eps0` and `adm_rescale()` refuses them.** `adm_eps0()`'s reference rule keeps
  individuals more than 99% of one non-focal ancestry; none of that arm's unadmixed individuals qualify,
  against 0.7% of the real cohort. Dividing two numbers selected to incomparable depths is not a rescale, and
  here it would compress the curve 23-fold -- the same error, repeated.
* Only the Indigenous-American estimand ships. The alternative, taking the whole AMR reference class as truth,
  has reliability that FALLS as true ancestry rises (0.83 at q = 0.05 to 0.40 at 0.39): a Latino founder's
  European segments count as American truth and are then correctly called European. That measures the
  reference class, not the caller.

# ADMIT 0.3.0

* `adm_threshold()` reports the PLATEAU, not only the optimum: every cut whose efficiency is within
  `plateau_tol` (default 0.02) of the best, which is the set the power criterion cannot tell apart. Reporting
  one cut invites the reader to treat it as special; in practice the efficiency curve is flat over a wide range
  -- on one real cohort it reaches 0.032 with no lower bound inside the evaluated range -- and hiding that hides the fact that the choice is
  still open. The object gains `plateau` and `conservative`.
* The plateau's upper end is named because the individuals it drops are the ones whose calls carry the most
  error, so it is the tightest cut power still treats as free. Whether to take it is a decision the package
  does not make: a cut that costs nothing in power can still cost sample size, and nothing in the efficiency
  calculation sees the calibration of the resulting scan.

# ADMIT 0.2.1

* `adm_eps0()` gains `background`: restrict the reference group to individuals near-pure in one named
  ancestry. `eps0` and the curve must be measured against the same background -- a reference group that is
  mostly European measures a European-background noise floor, and pairing it with an African-background curve
  calibrates one pairing with the other's error. Default behaviour is unchanged.

# ADMIT 0.2.0

* FLARE curves. `adm_curve()` gains `caller` (`"RFMix"`, `"FLARE"`) and `background`. The FLARE curves were
  measured on 16,000 genome-wide mosaics against a 3,201-sample 5-way reference matching the All of Us setup,
  with founders held out in five rotating folds, and were validated against the real 52,639-participant
  UK Biobank cohort run through the same FLARE setup: binned call variance agrees to within 15% across the
  whole range of global ancestry.
* A curve belongs to a PAIR of ancestries. South Asian ancestry on an African background is recovered almost
  perfectly (reliability 0.98 at q = 0.01), the same ancestry on a European background barely at all (0.09).
  Both are shipped for FLARE; pick the background your cohort has, or stratify.
* `adm_rescale()` warns beyond a twofold change in `eps0`. Rescaling slides the curve along q and cannot change
  the calling gain, so across a large gap it misstates reliability where q is large: rescaled to the South Asian
  European-background `eps0`, the African curve reports 0.94 at q = 0.6 where the measured value is 0.81.
* Shipped `eps0` values are now computed the way `adm_eps0()` computes them on real data -- the mean called focal
  ancestry among individuals more than 99% of one other ancestry -- so the two are on the same footing.
  A relatedness defect in the earlier simulation was fixed: 1000 Genomes trio children sat in the reference
  alongside their held-out parents, which inflated false calls severalfold.

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
