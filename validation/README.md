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
| `time_effects` | `survival::coxph(tt = )` | coefficients of the log-time form (veteran) |
| `adjusted_risk` | `riskRegression::ate`, weighted `survival::survfit` | g-formula and IPW risks |
| `turnbull` | `icenReg::ic_np`, `interval::icfit`, `interval::ictest` | NPMLE log-likelihood and survival, logrank test (Finkelstein scores, pclt) |
| `parsurv(Interval())` | `survival::survreg(interval2)`, `flexsurv::flexsurvreg` | coefficients, SEs, log-likelihood; exact times written as intervals give the `Event()` fit (also `directcif`) |
| `cox` / `parsurv` + `frailty()` | `frailtyEM::emfrail`, `parfm::parfm`, `survival::coxph(frailty())` | coefficients, frailty variance, LRT, log-likelihood (coxph's penalized fit reported only) |
| `directcif(estimation = "constrained", timing = )` | independent R log-likelihood; `stats::optim` (BFGS) of it | log-likelihood, coefficients, never below the nested model without `timing` |

Tolerances: 1e-6 for quantities computed by the same formula in both packages; 1e-3 (relative)
when the reference maximises the likelihood with a general-purpose optimiser; "info" when the two
packages use different valid estimators (reported, not judged).

`directcif` has no implementation in another package; it is validated in Part 3 against data
generated from the model. Part 3 also covers interval-censored Weibull regression (with the naive
fit to the midpoints for comparison) and shared gamma frailty (Cox and Weibull, coefficients and
log(theta), with the naive Cox fit that ignores the clusters), and the constrained direct model with
time ratios of the competing causes (coefficients, time ratios and the risk of a competing cause
at a fixed time, with the model without time ratios for comparison).

## Running

Install the reference packages once:

```r
install.packages(c("survival", "cmprsk", "riskRegression", "prodlim", "flexsurv", "flexsurvcure",
                   "survRM2", "timeROC", "rmarkdown", "icenReg", "interval", "frailtyEM", "parfm"))
# interval needs Icens from Bioconductor: install.packages("BiocManager"); BiocManager::install("Icens")
```

From this folder (or with the path to it):

```
Rscript run_validation.R                          # quick: every check once, small Monte Carlo (20-40 min)
Rscript run_validation.R full cores=6 run=v2      # full study for the paper (several hours)
Rscript run_validation.R full cores=6 run=v2 parts=1     # one part at a time, into the same run
Rscript run_validation.R run=v2 parts=report      # only remake the report from the finished stages
Rscript run_validation.R full sim_R=2000          # override any setting
```

The work is split into stages: accuracy has two (the established methods; the methods added
after 1.0.0), speed one per sample size and one for the number of covariates, the Monte Carlo
study one per sample size, stability one. Every finished stage is saved at once to
`validation_results/<run>/stages/`. A run can be stopped (Ctrl+C) and continued later with the same
command: finished stages are skipped (`redo=TRUE` runs them again), and only the stage that was
running is lost. A continued run must use the settings it was started with (they are kept in
`config.csv`; `force=TRUE` overrides). Every stage is seeded on its own, so the results do not
depend on how the run was split into sessions.

While a stage runs, a progress line is printed after every unit of work (a replication, a 2% share
of the Monte Carlo replications, a case of the stability checks):

```
14:05:12    [stage 4/12, run  28%] speed: n = 10000: 3/10 ( 30%), 27 min so far, about 63 min left
```

The share of the stage and the time left come from the units done so far; the share of the run
counts stages, which take very different times (the largest sample sizes take the longest).

Settings (`make_config()` in `R/00_common.R`): `sim_R` replications (quick 100, full 1000),
`sim_n` sample sizes (quick 300; full 200, 500, 1000), `speed_n` (quick 500, 2000, 10000; full
also 50000 and 200000), `speed_p`, `speed_reps`, `speed_max_seconds` (a method slower than this is
not run at larger n), `acc_reps`, `cores`, `seed`; `run`, `parts`, `redo`, `force`.

After every session the report is made again from all finished stages: the results file lists the
stages with the SurvX version, time and duration of each, says which stages are not finished yet,
warns when stages were run with different SurvX versions, and holds the log of every session.

Results go to `validation_results/<run>/` (without `run=`, a new run named by the date and time). The folder is not part of the package build.

## Quick checks of new functions

Short scripts, run from the repository root after `R CMD INSTALL .`, that print every result next
to the established package and write a text file and a PDF of the figures:

| Script | Functions | Compared with |
|---|---|---|
| `check_time_effects.R` | `time_effects()`, `adjusted_risk()` | `survival::coxph(tt = )` |
| `check_interval.R` | `Interval()`, `turnbull()`, `parsurv()` and `directcif()` with interval-censored data | `survival::survreg(type = "interval2")`, `interval::icfit`/`ictest`, `icenReg::ic_np`, `flexsurv::flexsurvreg` (when installed) |
| `check_frailty.R` | `cox()` and `parsurv()` with `frailty()` | `survival::coxph(frailty())`, `frailtyEM::emfrail`, `parfm::parfm` (when installed) |
| `check_assess.R` | `cindex()`, standard errors and paired differences of `brier_score()`/`auc_t()` and of the comparison tables, `calibration()` (O/E, ICI, E50, E90), `decision_curve()`, `anova()`, the new columns of `compare_models()`/`model_scan()` | `survival::concordance`, `coxph`, `survfit`, `finegray`, `survreg`, `anova.coxph`; a direct R computation of the C index over all pairs; `riskRegression::Score`, `pec::cindex`, `dcurves::dca` (when installed). About a minute; writes `check_assess.txt` with PASS/FAIL per comparison and the numbers of both packages for every disagreement |

## C++ memory check

`memcheck/` checks the whole C++ core with AddressSanitizer, UndefinedBehaviorSanitizer and
valgrind, without R; see `memcheck/README.md`.
