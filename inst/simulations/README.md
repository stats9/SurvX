# Simulation study: direct parametric CIF models versus Fine-Gray

**Question.** When do direct parametric models for the cumulative incidence function (CIF) predict
the absolute risk of an event better than the Fine-Gray model, and when do they not?

## Files

| File | Content |
| --- | --- |
| `sim_core.R` | Shared engine: scenarios, data generation, true CIF, model fitting, replications, summary. Sourced by the other scripts; not run by itself. |
| `sim_published_models.R` | **Study 1.** Published models only: direct PH / PO / GOR (Jeong & Fine 2006, 2007), the constrained model of Shi, Cheng & Jeong (2013), Fine-Gray, cause-specific Cox. |
| `sim_survx_models.R` | **Study 2.** Everything in study 1 plus the SurvX extensions: the weighted direct estimator, the constrained model with G chosen by AIC, spline baselines, and an automatic AIC choice. |
| `sim_figures.R` | Figures (PDF, 600-dpi TIFF, PNG) and tables (CSV, LaTeX) in BMC format, from the results of either or both studies. |

## How to run

Run everything from a terminal in this folder. In RStudio, set the working directory to this folder
and `source()` the scripts.

```sh
# 1. trial (a few minutes): checks that everything works on your machine
Rscript sim_published_models.R R=20 sizes=250
Rscript sim_figures.R

# 2. the real runs (hours; can be stopped and restarted, finished cells are kept)
Rscript sim_published_models.R R=1000 cores=8
Rscript sim_survx_models.R     R=1000 cores=8

# 3. figures and tables of both studies
Rscript sim_figures.R
```

**Optional settings:** `R=` replications per cell, `sizes=100,250,500,1000`, `scenarios=S1,S5`,
`cores=`, `n_test=`, `seed=`, `out_dir=`. The defaults are in the SETTINGS block at the top of
each script.

**Reproducibility:**
- Every replication has its own seed, so the results do not depend on the number of cores.
- The data sets of study 1 and study 2 are identical, so the models the two studies share give the same numbers.
- `results/setup_<study>.rds` stores the settings, the model list, the scenarios and the `sessionInfo()`.

**Run time:** study 2 fits 17 models per data set, and takes about three times as long as study 1.

## Design

### Data generation

**Covariates:** x1 ~ Bernoulli(0.5) and x2 ~ Uniform(−1, 1).

**Causes:** cause 1 is the cause of interest, and cause 2 competes with it.

**Censoring:** Uniform(0, c<sub>max</sub>), independent of everything else.

**Direct scenarios (S1–S6).** The true CIF is the Gompertz model of Jeong & Fine (2007, eq. 3.9):

F<sub>k</sub>(t; z) = 1 − {1 + α<sub>k</sub> e<sup>z′β<sub>k</sub></sup> u<sub>k</sub>(t)}<sup>−1/α<sub>k</sub></sup>,
u<sub>k</sub>(t) = τ<sub>k</sub>(e<sup>ρ<sub>k</sub>t</sup> − 1)/ρ<sub>k</sub>

- α = 0 gives proportional subdistribution hazards (PH), and α = 1 gives proportional odds (PO).
- Data are drawn in two steps:
  1. The cause is drawn from the plateaus F<sub>k</sub>(∞; z). With the remaining probability there is no event.
  2. The time is drawn from F<sub>k</sub>(t; z)/F<sub>k</sub>(∞; z).
- The code checks that F<sub>1</sub>(∞; z) + F<sub>2</sub>(∞; z) < 1 for every z.

### Why these scenarios

Under the PO and GOR models, the subdistribution hazard ratio of a covariate is not constant:
- It starts at e<sup>β</sup> and moves towards 1 as the cumulative incidence rises; the effect "wanes".
- Fine-Gray estimates one averaged ratio. Its predicted risks are therefore biased at some horizons, and this bias does not shrink as the sample grows.
- The departure from proportionality is small when the cumulative incidence stays low, because F/(1 − F) ≈ F.
- It is large when the incidence is high and the effects are strong.

`Fig2_shr_over_time` shows the ratio over time in every scenario.

The scenarios therefore:
- go from no departure (S1), through a small one (S2), to large ones (S3–S6, where α increases 1 → 2 → 4)
- include two scenarios that are neutral or unfavourable for the direct models (S7, S8)

The second point lets the study show both where the direct models gain and where they do not.

| | Truth | β cause 1; β cause 2 | Censoring | Horizons | Cause 1 / cause 2 / censored | Role |
| --- | --- | --- | --- | --- | --- | --- |
| S1 | direct PH (α = 0, 0) | (0.5, −0.4); (−0.3, 0.2) | U(0, 12) | 0.5, 1.5, 4 | 41 / 16 / 43 % | Fine-Gray correctly specified; only parametric efficiency can be gained |
| S2 | direct PO (α = 1, 0) | as S1 | U(0, 12) | 0.5, 1.5, 4 | 34 / 15 / 50 % | PO truth but low risks: mild non-proportionality |
| S3 | direct PO (α = 1, 1), strong effects | (1.5, −1.0); (−0.5, 0.3) | U(0, 12) | 0.5, 1.5, 4 | 55 / 9 / 36 % | clear non-proportionality |
| S4 | direct GOR (α = 2, 1) | (1.75, −1.25); (−0.5, 0.3) | U(0, 12) | 0.5, 1.5, 4 | 49 / 9 / 42 % | stronger waning |
| S5 | direct GOR (α = 4, 1) | (2.0, −1.5); (−0.5, 0.3) | U(0, 12) | 0.5, 1.5, 4 | 39 / 9 / 52 % | strongly waning effect |
| S6 | as S5, τ<sub>2</sub> = 0.1 | as S5 | U(0, 3) | 0.5, 1, 2 | 32 / 8 / 60 % | short follow-up, heavy censoring |
| S7 | Fine & Gray (1999) design, p = 0.4 | (0.5, −0.4); (−0.3, 0.3) | U(0, 5) | 0.3, 0.8, 2 | 40 / 40 / 20 % | Fine-Gray correct, Gompertz baseline wrong |
| S8 | exponential cause-specific hazards (0.5, 0.3) | (0.5, −0.4); (−0.3, 0.3) | U(0, 4) | 0.3, 0.8, 2 | 52 / 22 / 26 % | cause-specific Cox correct, everyone fails eventually |

**Further parameters.**
- S1–S5: ρ = (−0.5, −0.3).
- τ = (0.25, 0.08) in S1–S2 and (0.4, 0.05) in S3–S5.

**Event percentages** were computed from 200,000 simulated subjects. `sim_core.R` describes every
scenario in full, in the fields `role` and `why`. `sim_figures.R` writes these descriptions to
`tables/scenario_descriptions.txt`.

### Sample sizes and replications

- **Training sets:** n = 100, 250, 500 and 1000, with 1000 replications per cell.
- **Redraws:** a data set with fewer than 5 cause-1 or 3 cause-2 events is redrawn, and the number of redraws is reported.
- **Test set:** one per scenario, with 10,000 subjects and their true risks and uncensored outcomes.

### Performance measures

All measures are computed for the risk of cause 1 at three horizons:

| Measure | Definition |
| --- | --- |
| **RMSE** (primary) | √ mean over replications of the test-set mean of (predicted − true risk)². Censoring and the irreducible outcome noise do not affect it, so it is the most sensitive measure of prediction accuracy. |
| **% reduction against Fine-Gray** | 100 (1 − RMSE<sub>model</sub>/RMSE<sub>FG</sub>), computed on paired replications, with a 95 % Monte Carlo interval (delta method) |
| **win %** | the share of replications in which the model's error is below Fine-Gray's |
| **Bias** | mean signed error of the predicted risk |
| **Δ Brier** | difference in Brier score with the true outcomes of the test set (× 10⁴), with its Monte Carlo z. Differences between good models are small on this scale, because the Brier score is the outcome variance plus the MSE. |
| **Diagnostics** | convergence; the share of joint fits whose CIFs add up to more than 1; the model chosen by AIC; computing time |

**Model failures:** a fit that did not converge is still used for prediction, because dropping it
would favour the direct models. A fit that fails with an error is counted in `ok_pct`.

### What to expect

Pilot runs used the same C++ core with 30–60 replications per cell. They are only an indication;
the full study is the evidence.

- **S1, S2 (PH truth or low risks):** the direct models and Fine-Gray are within a few percent of each other.
- **Strong PO and GOR truths:** the RMSE reduction grows with α and with n:
  - S3, PO, n = 1000: about 30 % for the direct PO model
  - S5, α = 4: about 20–30 % at n = 300 and about 40–45 % at n = 1000 for the AIC-selected direct model
- **S7:** Fine-Gray is expected to be equal or better.
- **S8:** the joint direct model is expected to lose. The weighted and constrained estimators of study 2 are expected to close most of the gap.

## Figures and tables (`sim_figures.R`)

The figures follow BMC's figure guidelines: 85 mm or 170 mm wide, at most 225 mm high, vector PDF
plus 600-dpi TIFF with LZW compression. Text is 7–9 pt sans serif, and there are no titles inside
the figures. The Okabe-Ito colours come with distinct point shapes, so the figures also work in
grey scale.

| File | Content |
| --- | --- |
| Fig1_scenarios | true CIFs of both causes for four covariate profiles |
| Fig2_shr_over_time | subdistribution hazard ratio of x1 over time: why Fine-Gray is biased |
| Fig3_&lt;study&gt; | % RMSE reduction against Fine-Gray by n, middle horizon, with Monte Carlo intervals |
| FigS1_&lt;study&gt; | the same for every horizon (additional file) |
| Fig4_&lt;study&gt; | distribution over replications of the RMSE ratio to Fine-Gray, n = 500 |
| Fig5_&lt;study&gt; | bias of the predicted risk by horizon, n = 1000 |
| Fig6_published_vs_survx | automatic choices of the two studies side by side |
| Table1_scenarios | scenario definitions with the observed percentages |
| Table2_&lt;study&gt; | RMSE × 100 (MCSE) and change against Fine-Gray |
| TableS1_&lt;study&gt; | every measure for every cell (additional file) |

## References

- Beyersmann, J., Latouche, A., Buchholz, A. and Schumacher, M. (2009). Simulating competing risks data in survival analysis. *Statistics in Medicine*, 28, 956–971.
- Benichou, J. and Gail, M. H. (1990). Estimates of absolute cause-specific risk in cohort studies. *Biometrics*, 46, 813–826.
- Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the subdistribution of a competing risk. *JASA*, 94, 496–509.
- Jeong, J.-H. and Fine, J. (2006). Direct parametric inference for the cumulative incidence function. *JRSS C*, 55, 187–200.
- Jeong, J.-H. and Fine, J. P. (2007). Parametric regression on cumulative incidence function. *Biostatistics*, 8, 184–196.
- Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, 38, 2074–2102.
- Shi, H., Cheng, Y. and Jeong, J.-H. (2013). Constrained parametric model for simultaneous inference of two cumulative incidence functions. *Biometrical Journal*, 55, 82–96.
