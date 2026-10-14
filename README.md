# ADMIT

**AD**mixture-**M**apping **I**nclusion **T**hreshold

ADMIT chooses which individuals to include in an admixture-mapping scan. Instead of a fixed cutoff on
global ancestry (2%, 5%, self-identified population, or everyone), it derives the window from the
error rate of your local-ancestry caller and the ancestry distribution of your cohort, and reports
the power and effect-size attenuation of any cut you consider.

## Why an inclusion window

An individual contributes information to admixture mapping in proportion to the variance of their
local ancestry, 2q(1 − q), which vanishes as global ancestry q approaches 0 or 1. Local-ancestry
callers, however, assign a small amount of focal ancestry to people who carry none, at a background
rate ε₀. Near-unadmixed individuals therefore add calling noise but almost no signal, while a cutoff
that is too strict discards admixed individuals whose ancestry is called reliably. The power-maximizing
window lies between these extremes, at a few multiples of ε₀, and ADMIT locates it for your caller and
your cohort.

## Installation

```r
# install.packages("remotes")
remotes::install_github("AllyssaDrew/ADMIT")
```

ADMIT is written in base R (R ≥ 4.0) and has no compiled dependencies.

## Quick start

```r
library(ADMIT)

# 1. Background false-call rate of your caller, measured in your own cohort
eps <- adm_eps0(global_Q, local_dosage, focal = "AFR")

# 2. Inclusion window from a reliability curve measured against known truth
th <- adm_threshold(q_hat = global_Q$AFR, curve = "AFR", eps0 = eps)
th
plot(th)
```

`global_Q` holds each individual's global ancestry proportions (one column per ancestry), and
`local_dosage` is each individual's genome-wide mean focal local-ancestry dosage (0–2). Applied to
52,639 UK Biobank participants with RFMix calls, `adm_threshold()` returns:

```
Inclusion window  [measured curve]

  detection optimum    0.00382    efficiency 0.840   keeps 3,705
  estimation optimum   0.00568    efficiency 0.833   keeps 2,321

  recommended window   0.00382 to 0.00568
```

The **detection optimum** maximizes the non-centrality of the admixture-mapping test; the
**estimation optimum** minimizes the mean squared error of the effect estimate. Any cut between them
serves both purposes, and the efficiency curve is flat near its maximum, so the window can be
reported to one significant figure.

## Workflow

**Step 1. Measure ε₀.** `adm_eps0()` estimates the background rate from individuals who are
near-pure in a *non-focal* ancestry (default > 99%). Defining the reference group on the focal
estimate itself would select the individuals the caller got right and underestimate ε₀.

**Step 2. Choose a reliability curve.** ADMIT ships curves for African and South Asian focal ancestry,
measured on 16,000 simulated genomes called by RFMix v2 against a five-ancestry HGDP + 1000 Genomes
reference. Three routes are available:

| your setting | route |
|---|---|
| RFMix with a comparable reference | `curve = "AFR"` or `curve = "SAS"` |
| another caller | `adm_rescale(adm_curve("AFR"), eps0_to = eps)` |
| a new caller or reference panel | calibrate your own curve (below) |

Reliability depends on global ancestry through q/ε₀, so a shipped curve transfers to another caller
by rescaling with the two background rates. Between RFMix and FLARE, whose ε₀ differ 2.4-fold, the
rescaled curves agree to a median of 0.005; a FLARE user choosing the window from the rescaled RFMix
curve attains an efficiency of 0.858, against 0.859 from FLARE's own curve.

**Step 3. Compute the window.** `adm_threshold()` returns both optima, the efficiency and number of
individuals retained at every candidate cut, and the fraction of truly unadmixed individuals each cut
admits. `adm_efficiency()` evaluates any cut you specify.

## Calibrating a curve for your own caller

```r
split <- adm_holdout(ref_pop, frac = 0.25)            # founders the caller will not see

sim <- adm_simulate(hap, hap_pop, site, map,
                    q_dist = global_Q$AFR,            # your cohort's ancestry distribution
                    focal = "AFR", out_dir = "sim/")  # writes sim/sim_query.vcf

# run your caller on sim/sim_query.vcf against the reduced reference panel, then:
cur <- adm_score(sim, called = "sim/calls.msp.tsv", format = "msp")
adm_threshold(global_Q$AFR, curve = cur, eps0 = eps)
```

`adm_score()` reads RFMix and Gnomix `.msp.tsv` output, FLARE ancestry VCFs, or a pre-parsed matrix.
Simulate genome-wide: the window is set where focal tracts are rarest, and a single chromosome holds
too few of them to measure reliability there.

## Main functions

| function | purpose |
|---|---|
| `adm_eps0()` | background false-call rate of the caller, from your cohort |
| `adm_threshold()` | detection- and estimation-optimal inclusion window |
| `adm_efficiency()` | power efficiency and attenuation at any set of cuts |
| `adm_curve()` | shipped reliability curves (RFMix v2; AFR, SAS) |
| `adm_rescale()` | transfer a curve to another caller via ε₀ |
| `adm_holdout()`, `adm_simulate()`, `adm_score()` | calibrate a curve for a new caller |
| `adm_reliability()` | closed-form reliability as a function of q and ε₀ |

A worked example is in the vignette: `vignette("admit", package = "ADMIT")`.

## Citation

If you use ADMIT, please cite:

> Gan, Z. et al. Calibrated admixture mapping reveals ancestry-specific regulation of the human
> plasma proteome. *In preparation* (2026).

## Contact

Zeyu Gan (ganzeyu@unc.edu). Please report bugs and feature requests on the
[issue tracker](https://github.com/AllyssaDrew/ADMIT/issues).

## License

MIT
