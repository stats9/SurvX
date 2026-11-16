# SurvX validation suite

A reproducible study of SurvX in four parts. It answers two questions: is SurvX correct, stable
and fast enough to release, and what evidence goes into the software paper (JSS)? Every check is
recorded with its numbers, so a failure points to the exact function, data set and quantity.

| Part | What it does | Main output |
|---|---|---|
| 1. Accuracy | every SurvX estimate against the same quantity from an established package, on simulated designs (standard, tied times, delayed entry + strata, n = 5000) and on real data (lung, veteran, colon, pbc, mgus2) | `accuracy.csv`, `tables/accuracy_summary.*` |
| 2. Speed | time of one fit against n (500 to 200 000) and against the number of covariates, SurvX and reference on the same data | `speed.csv`, `tables/speed_by_n.*` |
| 3. Monte Carlo | data generated from each model with known parameters; bias, empirical SD, mean SE, SE/SD, coverage of 95% intervals, convergence; size and power of the tests | `simulation_raw.csv`, `simulation_summary.csv`, `tables/simulation.tex` |
| 4. Stability | difficult data (n = 15, 95% censoring, heavy ties, extreme scales, separation, one event, collinearity, no plateau) and invalid input; exact invariances (row order, rescaling, duplicated data, weights, seeds) | `stability.csv`, `tables/invariances.tex` |

Two files hold everything needed for review, small enough to send by message:

- `SurvX_validation_results.txt`: verdict, every summary table, every comparison that did not pass,
  the run log and the session information;
- `SurvX_validation_figures.pdf`: all figures (accuracy, coverage, SE/SD, speed).

The folder also keeps the complete CSV files of every check (for the archive, e.g. Zenodo), `tables/`
(CSV and LaTeX for the paper) and `validation_report.md` (and `.html` when rmarkdown and pandoc are
available).

## Comparisons (Part 1)

| SurvX | Reference | Quantities |
|---|---|---|
| `km` | `survival::survfit` | survival, Greenwood SE, median and CI |
| `logrank_test` | `survival::survdiff` | chi-square, rho = 0 and 1, stratified |
| `rmst` | `survival` (rmean), `survRM2::rmst2` | RMST and SE |
| `hazard_np` | `survival::survfit(ctype = 1)` | Nelson-Aalen |
| `aj`, `gray_test` | `cmprsk::cuminc`, `survival::survfit` (multi-state) | CIF, Gray test |
| `cox` | `survival::coxph`, `survfit.coxph`, `cox.zph` | coefficients, variance (model and robust), log-likelihood, predicted survival, residuals, PH test; Efron, Breslow, weights, strata, delayed entry, clusters, time-varying covariates |
| `fg` | `cmprsk::crr`, `survival::finegray` + `coxph` | coefficients, Fine-Gray variance, predicted CIF; delayed entry |
| `csc` | `survival::coxph`, `riskRegression::CSC` | coefficients, variance, predicted CIF, with strata |
| `parsurv` | `flexsurv::flexsurvreg`, `flexsurvspline`, `survival::survreg` | log-likelihood, predicted survival, AFT coefficients; 7 distributions, splines, delayed entry |
| `cure` | `flexsurvcure` | log-likelihood, predicted survival; mixture and non-mixture |
| `brier_score`, `auc_t` | `riskRegression::Score`, `timeROC` | Brier score, null model, AUC |

Tolerances: 1e-6 for quantities computed by the same formula in both packages; 1e-3 (relative)
when the reference maximises the likelihood with a general-purpose optimiser; "info" when the two
packages use different valid estimators (reported, not judged).

`directcif` has no implementation in another package; it is validated in Part 3 against data
generated from the model.

## Running

Install the reference packages once:

```r
install.packages(c("survival", "cmprsk", "riskRegression", "prodlim", "flexsurv", "flexsurvcure",
                   "survRM2", "timeROC", "rmarkdown"))
```

From this folder (or with the path to it):

```
Rscript run_validation.R                  # quick: every check once, small Monte Carlo (20-40 min)
Rscript run_validation.R full cores=6     # full study for the paper (several hours)
Rscript run_validation.R parts=1,4        # only accuracy and stability
Rscript run_validation.R full sim_R=2000  # override any setting
```

Settings (`make_config()` in `R/00_common.R`): `sim_R` replications (quick 100, full 1000),
`sim_n` sample sizes (quick 300; full 200, 500, 1000), `speed_n`, `speed_p`, `speed_reps`,
`speed_max_seconds` (a method slower than this is not run at larger n), `acc_reps`, `cores`, `seed`.

Results go to `validation_results/<date-time>/`. The folder is not part of the package build.
