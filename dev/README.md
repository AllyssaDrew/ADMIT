# dev/ — analysis scripts, not package tests

These build and validate the shipped curve against the UK Biobank data. They need `data.table` and
paths outside the package, so they must NOT live in `tests/`, where `R CMD check` would run them.

| script | what it does |
|---|---|
| `build_curve.R` | extracts the measured reliability curve into `R/sysdata.rda` |
| `rescale_test.R` | the RFMix-vs-FLARE test that validated `adm_rescale()` |
| `check_against_52k.R` | closed form vs the cohort; this is what showed the closed form is 5x too tight |
| `eps0_sensitivity.R` | how the optimum scales with eps0 |
| `validate_curve.R` | curve mode against the directly-computed windows |
| `check.sh` | `R CMD build` + `check` on a compute node with the MKL environment |
