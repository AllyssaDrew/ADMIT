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

  detection optimum    0.00198    efficiency 0.893   keeps 7,438
  estimation optimum   0.00241    efficiency 0.891   keeps 6,160

  recommended window   0.00198 to 0.00241   (anything inside is within reach of both)

  plateau              0.000688 to 0.00382   (efficiency within 0.02 of the optimum)
  conservative end     0.00382    efficiency 0.876   keeps 3,705
  Power cannot separate cuts inside the plateau. The conservative end costs at most
  0.02 of efficiency and drops 3,733 individuals whose calls are least reliable.
```

The **detection optimum** maximizes the non-centrality of the admixture-mapping test; the
**estimation optimum** minimizes the mean squared error of the effect estimate. Any cut between them
serves both purposes. The efficiency curve is flat near its maximum: the **plateau** is the range of
cuts within 0.02 of the optimum, and its upper end is the most conservative cut that power cannot
tell from the best one.

## Workflow

**Step 1. Measure ε₀.** `adm_eps0()` estimates the background rate from individuals who are
near-pure in a *non-focal* ancestry (default > 99%). Defining the reference group on the focal
estimate itself would select the individuals the caller got right and underestimate ε₀.

**Step 2. Choose a reliability curve.** ADMIT ships seven curves, each measured on simulated genomes of
known ancestry built genome-wide from individuals that the caller's reference panel does not contain:

| caller | reference panel | ancestry of interest (background) | simulated genomes |
|---|---|---|---|
| RFMix 2 | 446 individuals, five ancestries, HGDP + 1000 Genomes | African (European) | 5,489 |
| | | South Asian (European) | 5,563 |
| FLARE 0.6.0 | 3,201 individuals, five ancestries, HGDP + 1000 Genomes | African (European) | 6,500 |
| | | South Asian (European) | 6,486 |
| | | South Asian (African) | 2,997 |
| | | Indigenous American (European) | 2,967 |
| | | Indigenous American (African) | 904 |

A curve belongs to a caller, a reference panel and a pair of ancestries. Three routes are available:

| your setting | route |
|---|---|
| RFMix with a comparable reference | `curve = "AFR"` or `curve = "SAS"` |
| FLARE with a comparable reference | `adm_curve("AFR", caller = "FLARE")`, also `"SAS"` and `"AMR"`, with `background = "European"` or `"African"` |
| another caller with a similar ε₀ | `adm_rescale(adm_curve("AFR"), eps0_to = eps)` |
| a new caller or reference panel | calibrate your own curve (below) |

Rescaling slides a curve along q by the ratio of two background rates. It is an approximation. The
shipped RFMix and FLARE curves for African ancestry, whose ε₀ differ 2.3-fold, illustrate its reach:
on the UK Biobank distribution of global ancestry, and judged by FLARE's own curve, the window chosen
from the rescaled RFMix curve retains an efficiency of 0.908, against 0.914 from FLARE's own curve and
0.855 from the RFMix curve without rescaling. The rescaled curve still differs from the measured one
where global ancestry is lowest, so calibrate directly when the window matters.

**Step 3. Compute the window.** `adm_threshold()` returns both optima, the efficiency and number of
individuals retained at every candidate cut, and the fraction of truly unadmixed individuals each cut
admits. `adm_efficiency()` evaluates any cut you specify. Both accept `call_var`, the variance of each
individual's called dosage along the genome: above a global ancestry of 30% the curves assume more
variance than real calls show, and `call_var` replaces the assumed value by the observed one.

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
| `adm_curve()` | shipped reliability curves (RFMix 2 and FLARE 0.6.0; AFR, SAS, AMR; European or African background) |
| `adm_rescale()` | transfer a curve to another caller via ε₀ |
| `adm_holdout()`, `adm_simulate()`, `adm_score()` | calibrate a curve for a new caller |
| `adm_reliability()` | closed-form reliability as a function of q and ε₀ |

A worked example is in the vignette: `vignette("admit", package = "ADMIT")`.

## Citation

If you use ADMIT, please cite:

> Gan, Z. et al. Calibrated admixture mapping reveals ancestry differences in the plasma proteome.
> *In preparation* (2026).

## Contact

Zeyu Gan (ganzeyu@unc.edu). Please report bugs and feature requests on the
[issue tracker](https://github.com/AllyssaDrew/ADMIT/issues).

## License

MIT
