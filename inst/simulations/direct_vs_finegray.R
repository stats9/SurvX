#!/usr/bin/env Rscript
# When does a direct parametric model for the cumulative incidence function predict better than
# the Fine-Gray model? Six data-generating scenarios, three sample sizes.
#
#   Rscript direct_vs_finegray.R [replications = 200] [sample sizes = 100,300,1000]
#
# Design
#   Training data of size n with uniform censoring; every model is fitted to the training data and
#   predicts the risk of cause 1 at three horizons for 20000 independent test subjects.
#   Two measures on the test subjects, both free of censoring because the data are simulated:
#     rmse   root mean squared difference between the predicted and the TRUE risk F_1(t; z)
#     brier  mean squared difference between the predicted risk and the event indicator
#   The Brier score equals (irreducible outcome variance) + (mean squared error of the risk), so
#   differences between good models are small on the Brier scale; they are reported as paired
#   differences to Fine-Gray with a Monte Carlo z statistic.
#   "direct AIC" is the direct model (PH, PO or GOR link) with the smallest AIC in each replication,
#   i.e. what an analyst would do without knowing the truth.
#
# Scenarios
#   A  truth: direct Gompertz model, proportional subdistribution hazards (PH). Fine-Gray is
#      correctly specified too; the parametric model can only gain efficiency.
#   B  truth: direct Gompertz model, proportional odds (alpha = 1), strong covariate effects.
#      The subdistribution hazards are not proportional: Fine-Gray is misspecified.
#   C  as A with heavy early censoring and a horizon near the end of follow-up.
#   D  truth: exponential cause-specific hazards. Every model here is misspecified.
#   E  truth: generalized odds rate model with alpha = 4 (strongly non-proportional).
#   F  truth: the simulation design of Fine & Gray (1999, Section 5):
#      F_1(t; z) = 1 - {1 - p (1 - e^-t)}^exp(z'beta_1), cause 2 exponential given that it occurs.
#      Proportional subdistribution hazards with a baseline that is NOT Gompertz: Fine-Gray is
#      correctly specified, the direct Gompertz models are not. A neutral scenario.
#
# Models (all with both covariates; direct models with a Gompertz baseline)
#   joint PH / PO / GOR   all causes by the full likelihood of Jeong & Fine (2006, 2007)
#   joint AIC             the joint model with the smallest AIC
#   weighted PH / PO / GOR  cause 1 alone, censoring weights as in Fine & Gray (1999); directcif(estimation = "weighted")
#   weighted AIC          the weighted model with the smallest pseudo-AIC
#   constrained PH / PO   Shi, Cheng & Jeong (2013): asymptotes add up to one; G chosen by AIC among
#                         the modified logistic, Weibull, log-logistic and log-normal forms
#   direct auto           smallest AIC among the joint and the constrained models (both are full
#                         likelihoods of the same data, so their AIC values are comparable)
#   Fine-Gray, cause-specific Cox
# A fit that did not converge is still used for prediction (dropping it would favour the direct
# models); it is only excluded from the AIC choices when a converged candidate exists.
suppressPackageStartupMessages(library(SurvX))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 200L
sizes <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(100L, 300L, 1000L)
set.seed(20261114)

scenarios <- list(
  A_PH = list(kind = "direct", rho = c(-0.5, -0.3), tau = c(0.25, 0.08), alpha = c(0, 0),
              beta = list(c(0.5, -0.4), c(-0.3, 0.2)), cmax = 12, times = c(0.5, 1.5, 4)),
  B_PO = list(kind = "direct", rho = c(-0.5, -0.3), tau = c(0.4, 0.05), alpha = c(1, 1),
              beta = list(c(1.5, -1.0), c(-0.5, 0.3)), cmax = 12, times = c(0.5, 1.5, 4)),
  C_PH_heavy_censoring = list(kind = "direct", rho = c(-0.5, -0.3), tau = c(0.25, 0.08), alpha = c(0, 0),
              beta = list(c(0.5, -0.4), c(-0.3, 0.2)), cmax = 2, times = c(0.5, 1.0, 1.8)),
  D_cause_specific_hazards = list(kind = "csh", h = c(0.5, 0.3), beta = list(c(0.5, -0.4), c(-0.3, 0.3)),
              cmax = 4, times = c(0.3, 0.8, 2.0)),
  E_GOR_alpha4 = list(kind = "direct", rho = c(-0.5, -0.3), tau = c(0.4, 0.05), alpha = c(4, 1),
              beta = list(c(2.0, -1.5), c(-0.5, 0.3)), cmax = 12, times = c(0.5, 1.5, 4)),
  F_FineGray_design = list(kind = "fg", p = 0.4, beta = list(c(0.5, -0.4), c(-0.3, 0.3)), cmax = 5, times = c(0.3, 0.8, 2.0)))

# F(t; z) = 1 - {1 + alpha e^eta u(t)}^(-1/alpha), u(t) = tau (e^(rho t) - 1) / rho   (Jeong & Fine 2007, eq. 3.9)
Fgor <- function(t, eta, rho, tau, alpha) {
  u <- ifelse(is.finite(t), tau * expm1(rho * t) / rho, -tau / rho)
  eu <- exp(eta) * u
  if (alpha == 0) -expm1(-eu) else 1 - (1 + alpha * eu)^(-1 / alpha)
}
Finv <- function(Fv, eta, rho, tau, alpha) {
  eu <- if (alpha == 0) -log1p(-Fv) else ((1 - Fv)^(-alpha) - 1) / alpha
  log1p(rho * eu / exp(eta) / tau) / rho
}
gen <- function(n, S) {
  X <- cbind(x1 = rbinom(n, 1, 0.5), x2 = runif(n, -1, 1))
  if (S$kind == "direct") {
    eta <- lapply(S$beta, function(b) drop(X %*% b))
    pl <- lapply(1:2, function(k) Fgor(rep(Inf, n), eta[[k]], S$rho[k], S$tau[k], S$alpha[k]))
    stopifnot(max(pl[[1]] + pl[[2]]) < 1)
    u <- runif(n); cause <- ifelse(u < pl[[1]], 1L, ifelse(u < pl[[1]] + pl[[2]], 2L, 0L)); Tt <- rep(Inf, n)
    for (k in 1:2) { i <- cause == k; Tt[i] <- Finv(runif(sum(i)) * pl[[k]][i], eta[[k]][i], S$rho[k], S$tau[k], S$alpha[k]) }
  } else if (S$kind == "fg") {
    e <- exp(drop(X %*% S$beta[[1]])); P1 <- 1 - (1 - S$p)^e
    cause <- ifelse(runif(n) < P1, 1L, 2L); Fv <- runif(n) * P1
    t1 <- -log(1 - (1 - (1 - Fv)^(1 / e)) / S$p); t2 <- rexp(n, exp(drop(X %*% S$beta[[2]])))
    Tt <- ifelse(cause == 1L, t1, t2)
  } else {
    t1 <- rexp(n, S$h[1] * exp(drop(X %*% S$beta[[1]]))); t2 <- rexp(n, S$h[2] * exp(drop(X %*% S$beta[[2]])))
    Tt <- pmin(t1, t2); cause <- ifelse(t1 < t2, 1L, 2L)
  }
  list(X = X, T = Tt, cause = cause)
}
truth <- function(X, t, S) {
  if (S$kind == "direct") return(Fgor(rep(t, nrow(X)), drop(X %*% S$beta[[1]]), S$rho[1], S$tau[1], S$alpha[1]))
  if (S$kind == "fg") return(1 - (1 - S$p * (1 - exp(-t)))^exp(drop(X %*% S$beta[[1]])))
  h1 <- S$h[1] * exp(drop(X %*% S$beta[[1]])); h2 <- S$h[2] * exp(drop(X %*% S$beta[[2]]))
  h1 / (h1 + h2) * -expm1(-(h1 + h2) * t)
}

models <- c("joint PH", "joint PO", "joint GOR", "joint AIC", "weighted PH", "weighted PO", "weighted GOR", "weighted AIC",
            "constrained PH", "constrained PO", "direct auto", "Fine-Gray", "cause-specific Cox")
links <- c("PH", "PO", "GOR"); conds <- c("logistic", "weibull", "loglogistic", "lognormal")
fml <- Event(time, status) ~ x1 + x2
aic_of <- function(f) -2 * f@loglik + 2 * length(f@coefficients)
best_of <- function(fits) {                 # smallest AIC, converged fits first
  ok <- vapply(fits, function(f) isTRUE(f@converged), logical(1))
  if (any(ok)) fits <- fits[ok]
  fits[[which.min(vapply(fits, aic_of, numeric(1)))]]
}
quiet <- function(expr) tryCatch(expr, error = function(e) NULL)

results <- NULL
for (sn in names(scenarios)) for (n in sizes) {
  S <- scenarios[[sn]]; nt <- length(S$times)
  te <- gen(20000L, S); newd <- as.data.frame(te$X)
  tr <- sapply(S$times, function(s) truth(te$X, s, S))
  Y <- sapply(S$times, function(s) as.numeric(te$T <= s & te$cause == 1L))
  mse <- br <- array(NA_real_, c(R, nt, length(models)), dimnames = list(NULL, NULL, models))
  conv <- stats::setNames(numeric(length(models)), models); over1 <- conv
  chosen <- list("joint AIC" = table(factor(character(0), links)), "weighted AIC" = table(factor(character(0), links)),
                 "direct auto" = table(factor(character(0), c("joint", "constrained"))))
  t0 <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    g <- gen(n, S); cens <- runif(n, 0, S$cmax)
    d <- data.frame(time = pmin(g$T, cens), status = ifelse(g$T <= cens, g$cause, 0L), g$X)
    if (sum(d$status == 1) < 8 || sum(d$status == 2) < 4) next
    d$time <- pmax(d$time, 1e-8)
    fj <- lapply(links, function(lk) quiet(directcif(fml, d, link = lk)))
    fw <- lapply(links, function(lk) quiet(directcif(fml, d, link = lk, estimation = "weighted", cause = 1)))
    fc <- lapply(c("PH", "PO"), function(lk) {
      cand <- Filter(Negate(is.null), lapply(conds, function(cd) quiet(directcif(fml, d, link = lk, estimation = "constrained", cause = 1, cond = cd))))
      if (length(cand)) best_of(cand) else NULL })
    ffg <- quiet(fg(fml, d)); fcs <- quiet(csc(fml, d))
    if (any(vapply(c(fj, fw, fc), is.null, logical(1))) || is.null(ffg) || is.null(fcs)) next
    names(fj) <- names(fw) <- links
    bj <- best_of(fj); bw <- best_of(fw); ba <- best_of(c(fj, fc))
    fits <- c(stats::setNames(fj, paste("joint", links)), list("joint AIC" = bj), stats::setNames(fw, paste("weighted", links)),
              list("weighted AIC" = bw, "constrained PH" = fc[[1]], "constrained PO" = fc[[2]], "direct auto" = ba))
    chosen[["joint AIC"]][bj@link[1]] <- chosen[["joint AIC"]][bj@link[1]] + 1L
    chosen[["weighted AIC"]][bw@link[1]] <- chosen[["weighted AIC"]][bw@link[1]] + 1L
    chosen[["direct auto"]][ba@estimation] <- chosen[["direct auto"]][ba@estimation] + 1L
    P <- lapply(fits, function(f) predict(f, newd, times = S$times, cause = 1))
    P[["Fine-Gray"]] <- predict(ffg, newd, times = S$times)
    P[["cause-specific Cox"]] <- predict(fcs, newd, times = S$times, cause = 1)
    for (m in names(fits)) { conv[m] <- conv[m] + isTRUE(fits[[m]]@converged); over1[m] <- over1[m] + (fits[[m]]@diagnostics$n_exceed_range > 0) }
    for (m in models) { mse[r, , m] <- colMeans((P[[m]] - tr)^2); br[r, , m] <- colMeans((P[[m]] - Y)^2) }
  }
  keep <- !is.na(mse[, 1, "Fine-Gray"]); nrep <- sum(keep)
  direct <- setdiff(models, c("Fine-Gray", "cause-specific Cox"))
  for (m in models) for (j in seq_len(nt)) {
    dlt <- br[keep, j, m] - br[keep, j, "Fine-Gray"]
    results <- rbind(results, data.frame(
      scenario = sn, n = n, time = S$times[j], model = m, replications = nrep,
      rmse_vs_true_risk = sqrt(mean(mse[keep, j, m])), median_rmse = median(sqrt(mse[keep, j, m])), brier = mean(br[keep, j, m]),
      brier_minus_FG_x1e4 = 1e4 * mean(dlt), z = if (m == "Fine-Gray") NA else mean(dlt) / (sd(dlt) / sqrt(nrep)),
      pct_replications_better_than_FG = if (m == "Fine-Gray") NA else 100 * mean(dlt < 0),
      converged_pct = if (m %in% direct) 100 * conv[[m]] / nrep else NA,
      sum_of_CIFs_above_1_pct = if (m %in% direct) 100 * over1[[m]] / nrep else NA,
      choice = if (m %in% names(chosen)) paste(names(chosen[[m]]), as.integer(chosen[[m]]), collapse = " ") else NA))
  }
  cat("\n====", sn, " n =", n, " replications =", nrep, " seconds =", round(proc.time()[["elapsed"]] - t0), "====\n")
  print(results[results$scenario == sn & results$n == n, c("time", "model", "rmse_vs_true_risk", "brier_minus_FG_x1e4", "z", "converged_pct")], row.names = FALSE, digits = 3)
  write.csv(results, "direct_vs_finegray_results.csv", row.names = FALSE)     # saved after every cell
}
cat("\nwritten: direct_vs_finegray_results.csv\n",
    "rmse_vs_true_risk: smaller is better (median_rmse: the same per replication, robust to a few bad fits).\n",
    "brier_minus_FG_x1e4 < 0 and z < -2: the model beats Fine-Gray.\n", sep = "")
