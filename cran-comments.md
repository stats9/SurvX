## Submission

This is the first submission of SurvX to CRAN.

SurvX provides estimators, tests and regression models for survival and competing risks data
(Kaplan-Meier, Aalen-Johansen, log-rank and Gray tests, Cox, cause-specific Cox, Fine-Gray,
parametric, cure and direct parametric models for the cumulative incidence function) in a C++
core without external C++ dependencies, with an S7 interface.

## Test environments

* GitHub Actions: Ubuntu (R release and R devel), macOS (R release), Windows (R release)
* GitHub Actions: `R CMD check --as-cran` with `_R_CHECK_CRAN_INCOMING_` on Ubuntu, R release and
  R devel

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  Words flagged as possibly misspelled in DESCRIPTION, if any, are names of methods or authors
  (Aalen, Johansen, Kaplan, Meier, Nelson, Gompertz, Royston, Parmar, Fine, Gray, Jeong) and are
  spelled correctly.

## Notes for the reviewers

* The package uses `.C()` registered routines (see `src/init.cpp`); all inputs are validated in R
  and again in C++.
* Every function of the package was compared with the established reference implementations
  (survival, cmprsk, riskRegression, flexsurv, flexsurvcure, survRM2); the validation suite is in
  the GitHub repository (directory `validation`) and is excluded from the package.
* The data set `fertility` is also distributed with the CRAN package cmpp by the same author;
  its source publication is cited in `?fertility`.
