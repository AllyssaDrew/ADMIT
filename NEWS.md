# ADMIT 0.3.0

* `adm_threshold()` reports the PLATEAU, not only the optimum: every cut whose efficiency is within
  `plateau_tol` (default 0.02) of the best, which is the set the power criterion cannot tell apart. Reporting
  one cut invites the reader to treat it as special; in practice the efficiency curve is flat over a wide range
  -- for one real cohort, 6e-5 to 0.06, a thousandfold -- and hiding that hides the fact that the choice is
  still open. The object gains `plateau` and `conservative`.
* The plateau's upper end is named because the individuals it drops are the ones whose calls carry the most
  error, so it is the tightest cut power still treats as free. Whether to take it is a decision the package
  does not make: a cut that costs nothing in power can still cost sample size, and nothing in the efficiency
  calculation sees the calibration of the resulting scan.

