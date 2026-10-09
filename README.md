# ADMIT — ADMixture-mapping Inclusion Threshold

Chooses the global-ancestry window for an admixture-mapping cohort from the local-ancestry caller's
own error rate, instead of by convention.

```r
eps <- adm_eps0(global_Q, local_dosage, focal = "AFR")   # your caller's floor, from your own cohort
th  <- adm_threshold(cohort$AFR, curve = "AFR", eps0 = eps)
plot(th)
```

## Why

An admixture-mapping scan's non-centrality parameter is `(sum a)^2 / sum w` — the signal each
individual contributes over the variance each contributes. An individual raises it only when their
call is more reliable than about half the pooled reliability. Reliability is
`2q(1-q) / (2q(1-q) + eps0)`: the numerator vanishes as an individual approaches unadmixed, the
caller's noise does not. So a threshold has to exist, and where it goes is set by `eps0` — which you
can measure in any cohort, because it is what the caller reports in people who carry none of the
focal ancestry.

## Status and honest limits

- **`eps0` estimation is direct and needs no simulation.** Define the reference group by being
  near-pure in some *other* ancestry; defining it on the focal estimate is circular.
- **The one-parameter closed form gets the scaling right and the level wrong.** Against known truth
  it places the optimum ~5x too tight and overstates efficiency by ~0.18. `adm_threshold()` warns if
  you use it without a curve.
- **Rescaling a shipped curve to another caller works, and is tested.** RFMix 2 and FLARE on the
  same 16,000 mosaics differ 2.4-fold in `eps0`; against the scaled axis `q/eps0` their measured
  reliability curves agree to a median of 0.005. On a real 52,639-individual cohort a FLARE user
  attains efficiency 0.859 from FLARE's own curve, **0.858** from the rescaled RFMix curve, 0.843
  from the unrescaled one. Use `adm_rescale()`.
- **The choice is forgiving.** The efficiency curve is a plateau: being twofold off on the cut costs
  under 0.02. Report the window to one significant figure.
- Validated on one caller pair, one focal ancestry, one reference panel.
- **Calibrating your own curve is wrapped** (`adm_holdout` -> `adm_simulate` -> your caller ->
  `adm_score`). The package generates the query and scores the answer; it never runs your caller,
  because that is the step whose command lines and configs differ between tools and releases.
  `adm_score()` reads the `.msp` layout that RFMix and Gnomix share, or takes a decoded matrix.
- **Calibrate genome-wide.** The low-`q` end sets the threshold and needs the most data. On chr22
  alone, reliability matches the genome-wide answer to 0.016 above `q_hat` 0.025 but collapses
  below it (0.13 against 0.84). `adm_score()` warns when the simulation is too thin.

## Install

```r
remotes::install_github("AllyssaDrew/ADMIT")
```
