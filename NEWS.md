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

