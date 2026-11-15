# =====================================================================================================
#  sim_core.R -- shared engine of the simulation study
#  "When do direct parametric models for the cumulative incidence function predict better than
#   the Fine-Gray model?"
#
#  This file is sourced by sim_published_models.R and sim_survx_models.R; it is not run by itself.
#  It holds
#    1. the data-generating scenarios  (scenarios(), with a full description of every parameter)
#    2. data generation and the true cumulative incidence  (gen_data(), true_cif())
#    3. model fitting and prediction for a list of model specifications  (fit_models())
#    4. one replication, one cell (scenario x sample size), and the whole study with resumption
#       and parallel computing  (one_replication(), run_cell(), run_study())
#    5. the summary with Monte Carlo standard errors  (summarise_study())
#
#  Requirements: R >= 4.1, SurvX >= 0.5.1; the package 'parallel' (part of base R).
# =====================================================================================================

suppressPackageStartupMessages(library(SurvX))
if (utils::packageVersion("SurvX") < "0.5.1")
  stop("the simulation scripts need SurvX 0.5.1 or later; restart R and install the current version")

# -----------------------------------------------------------------------------------------------------
# 1. Scenarios
# -----------------------------------------------------------------------------------------------------
# Two covariates in every scenario: x1 ~ Bernoulli(0.5) (e.g. treatment or a binary risk factor) and
# x2 ~ Uniform(-1, 1) (a centred continuous marker). Cause 1 is the cause of interest, cause 2 the
# competing cause. Time is in years. Censoring is Uniform(0, cmax), independent of everything else.
#
# Direct (Gompertz / generalised odds rate) truth, Jeong & Fine (2007, eq. 3.9):
#   F_k(t; z) = 1 - {1 + alpha_k exp(z'beta_k) u_k(t)}^(-1/alpha_k),  u_k(t) = tau_k {exp(rho_k t) - 1} / rho_k
#   alpha_k = 0 means 1 - exp{-exp(z'beta_k) u_k(t)} (proportional subdistribution hazards, PH);
#   alpha_k = 1 is the proportional odds (PO) model. rho_k < 0 gives a plateau F_k(inf; z) < 1.
#   Data: the cause is drawn from the plateaus (cause 1 with probability F_1(inf; z), cause 2 with
#   F_2(inf; z), no event otherwise, which requires F_1(inf; z) + F_2(inf; z) < 1 for every z; this
#   is checked), then the event time from F_k(t; z) / F_k(inf; z).
#
# Why these scenarios. Under the PO and GOR models the subdistribution hazard ratio of a covariate is
# not constant: it starts at exp(beta) and moves towards 1 as the cumulative incidence rises, i.e.
# the effect "wanes" on the hazard scale. A Fine-Gray model then averages a time-varying ratio and its
# predicted risks are biased at some horizons; the bias does not shrink with the sample size. The
# departure from proportionality is small when the cumulative incidence stays low (the odds F/(1-F)
# is then close to F) and large when it is high and the covariate effects are strong. The scenarios
# therefore move from "no departure" (S1) through "small" (S2) to "large" (S3-S6), and include two
# scenarios in which the Fine-Gray model is correctly specified or in which no model is (S7, S8), so
# that the results show both where the direct models gain and where they do not.
scenarios <- function() list(
  S1 = list(id = "S1", label = "S1: PH, moderate effects", kind = "direct",
            rho = c(-0.5, -0.3), tau = c(0.25, 0.08), alpha = c(0, 0),
            beta = list(c(0.5, -0.4), c(-0.3, 0.2)), cmax = 12, times = c(0.5, 1.5, 4),
            role = "Fine-Gray correctly specified (proportional subdistribution hazards); the direct PH model is the true model. Any gain is parametric efficiency only.",
            why = "Plateau of cause 1 0.39 at x = 0 (0.28 to 0.71 over the covariate range), moderate effects (sHR 1.65 for x1, 0.67 per unit of x2); about 41% cause-1 events, 16% cause-2 events, 43% censored."),
  S2 = list(id = "S2", label = "S2: PO, moderate effects", kind = "direct",
            rho = c(-0.5, -0.3), tau = c(0.25, 0.08), alpha = c(1, 0),
            beta = list(c(0.5, -0.4), c(-0.3, 0.2)), cmax = 12, times = c(0.5, 1.5, 4),
            role = "Proportional odds truth with moderate effects and plateaus of cause 1 between 0.25 and 0.55: the subdistribution hazards are only mildly non-proportional.",
            why = "Same baseline and effects as S1, only the link changes; shows that a PO truth alone does not make Fine-Gray poor when risks are low. About 34% cause-1 events, 15% cause-2 events, 50% censored."),
  S3 = list(id = "S3", label = "S3: PO, strong effects", kind = "direct",
            rho = c(-0.5, -0.3), tau = c(0.4, 0.05), alpha = c(1, 1),
            beta = list(c(1.5, -1.0), c(-0.5, 0.3)), cmax = 12, times = c(0.5, 1.5, 4),
            role = "Proportional odds truth with strong effects (OR 4.5 for x1, 0.37 per unit of x2) and plateaus of cause 1 between 0.23 and 0.91: clearly non-proportional subdistribution hazards.",
            why = "Plateau of cause 1 0.44 at x = 0; about 55% cause-1 events, 9% cause-2 events, 36% censored."),
  S4 = list(id = "S4", label = "S4: GOR (alpha = 2), strong effects", kind = "direct",
            rho = c(-0.5, -0.3), tau = c(0.4, 0.05), alpha = c(2, 1),
            beta = list(c(1.75, -1.25), c(-0.5, 0.3)), cmax = 12, times = c(0.5, 1.5, 4),
            role = "Generalised odds rate truth between S3 and S5: stronger waning of the subdistribution hazard ratio than under PO.",
            why = "Intermediate point of the gradient alpha = 1, 2, 4 (S3, S4, S5), with effects between those of S3 and S5. Plateaus of cause 1 between 0.17 and 0.83; about 49% cause-1 events, 9% cause-2 events, 42% censored."),
  S5 = list(id = "S5", label = "S5: GOR (alpha = 4), strong effects", kind = "direct",
            rho = c(-0.5, -0.3), tau = c(0.4, 0.05), alpha = c(4, 1),
            beta = list(c(2.0, -1.5), c(-0.5, 0.3)), cmax = 12, times = c(0.5, 1.5, 4),
            role = "Strongly non-proportional subdistribution hazards: the effect of x1 is large early and fades later (a 'waning' effect).",
            why = "Plateaus of cause 1 between 0.13 and 0.69; about 39% cause-1 events, 9% cause-2 events, 52% censored. In pilot runs (30 replications) the direct GOR model had about 45% lower RMSE than Fine-Gray at n = 1000."),
  S6 = list(id = "S6", label = "S6: GOR (alpha = 4), heavy censoring", kind = "direct",
            rho = c(-0.5, -0.3), tau = c(0.4, 0.1), alpha = c(4, 1),
            beta = list(c(2.0, -1.5), c(-0.5, 0.3)), cmax = 3, times = c(0.5, 1, 2),
            role = "As S5 with censoring Uniform(0, 3) (short follow-up; horizons up to 2/3 of the maximum follow-up) and a competing cause twice as frequent (tau_2 = 0.1).",
            why = "Tests whether the advantage survives when most subjects are censored before the plateau is reached. About 32% cause-1 events, 8% cause-2 events, 60% censored."),
  S7 = list(id = "S7", label = "S7: Fine-Gray design", kind = "fg", p = 0.4,
            beta = list(c(0.5, -0.4), c(-0.3, 0.3)), cmax = 5, times = c(0.3, 0.8, 2.0),
            role = "Fine-Gray correctly specified with a baseline that is NOT Gompertz: every direct Gompertz model is misspecified. Neutral/unfavourable scenario for the direct models.",
            why = "The simulation design of Fine & Gray (1999, Section 5): F_1(t; z) = 1 - {1 - p(1 - exp(-t))}^exp(z'beta_1), p = 0.4; cause 2 exponential with rate exp(z'beta_2) given that it occurs. About 40% cause-1 events, 40% cause-2 events, 20% censored."),
  S8 = list(id = "S8", label = "S8: cause-specific hazards", kind = "csh", h = c(0.5, 0.3),
            beta = list(c(0.5, -0.4), c(-0.3, 0.3)), cmax = 4, times = c(0.3, 0.8, 2.0),
            role = "Exponential proportional cause-specific hazards: the cause-specific Cox model is correctly specified; Fine-Gray and every direct model are misspecified, and every subject eventually fails.",
            why = "A common data-generating mechanism in methodological studies (e.g. Beyersmann et al. 2009); the joint direct model assumes that some subjects never fail. About 52% cause-1 events, 22% cause-2 events, 26% censored.")
)

# -----------------------------------------------------------------------------------------------------
# 2. Data generation and truth
# -----------------------------------------------------------------------------------------------------
F_gor <- function(t, eta, rho, tau, alpha) {
  u <- ifelse(is.finite(t), tau * expm1(rho * t) / rho, -tau / rho)
  eu <- exp(eta) * u
  if (alpha == 0) -expm1(-eu) else 1 - (1 + alpha * eu)^(-1 / alpha)
}
F_gor_inverse <- function(Fv, eta, rho, tau, alpha) {          # time at which F_gor equals Fv
  eu <- if (alpha == 0) -log1p(-Fv) else ((1 - Fv)^(-alpha) - 1) / alpha
  log1p(rho * eu / exp(eta) / tau) / rho
}
gen_covariates <- function(n) cbind(x1 = stats::rbinom(n, 1L, 0.5), x2 = stats::runif(n, -1, 1))

gen_data <- function(n, S) {
  X <- gen_covariates(n)
  if (S$kind == "direct") {
    eta <- lapply(S$beta, function(b) drop(X %*% b))
    plat <- lapply(1:2, function(k) F_gor(rep(Inf, n), eta[[k]], S$rho[k], S$tau[k], S$alpha[k]))
    if (max(plat[[1]] + plat[[2]]) >= 1) stop("scenario ", S$id, ": plateaus add up to 1 or more")
    u <- stats::runif(n)
    cause <- ifelse(u < plat[[1]], 1L, ifelse(u < plat[[1]] + plat[[2]], 2L, 0L))
    Tt <- rep(Inf, n)
    for (k in 1:2) {
      i <- cause == k
      Tt[i] <- F_gor_inverse(stats::runif(sum(i)) * plat[[k]][i], eta[[k]][i], S$rho[k], S$tau[k], S$alpha[k])
    }
  } else if (S$kind == "fg") {
    e <- exp(drop(X %*% S$beta[[1]])); P1 <- 1 - (1 - S$p)^e
    cause <- ifelse(stats::runif(n) < P1, 1L, 2L)
    Fv <- stats::runif(n) * P1
    t1 <- -log(1 - (1 - (1 - Fv)^(1 / e)) / S$p)
    t2 <- stats::rexp(n, exp(drop(X %*% S$beta[[2]])))
    Tt <- ifelse(cause == 1L, t1, t2)
  } else {
    t1 <- stats::rexp(n, S$h[1] * exp(drop(X %*% S$beta[[1]])))
    t2 <- stats::rexp(n, S$h[2] * exp(drop(X %*% S$beta[[2]])))
    Tt <- pmin(t1, t2); cause <- ifelse(t1 < t2, 1L, 2L)
  }
  cens <- if (is.finite(S$cmax)) stats::runif(n, 0, S$cmax) else rep(Inf, n)
  data.frame(time = pmax(pmin(Tt, cens), 1e-8), status = ifelse(Tt <= cens, cause, 0L), X)
}

# true cumulative incidence of cause k at time t for covariate matrix X (n x 2); t may be Inf
true_cif <- function(X, t, S, k = 1L) {
  eta <- drop(X %*% S$beta[[k]])
  if (S$kind == "direct") return(F_gor(rep(t, nrow(X)), eta, S$rho[k], S$tau[k], S$alpha[k]))
  if (S$kind == "fg") {
    if (k == 1L) return(1 - (1 - S$p * (1 - exp(-t)))^exp(eta))
    P1 <- 1 - (1 - S$p)^exp(drop(X %*% S$beta[[1]]))
    return((1 - P1) * (1 - exp(-exp(eta) * t)))
  }
  h1 <- S$h[1] * exp(drop(X %*% S$beta[[1]])); h2 <- S$h[2] * exp(drop(X %*% S$beta[[2]]))
  hk <- if (k == 1L) h1 else h2
  hk / (h1 + h2) * -expm1(-(h1 + h2) * t)
}

# -----------------------------------------------------------------------------------------------------
# 3. Models
# -----------------------------------------------------------------------------------------------------
# A model specification is a list with
#   label   name shown in tables and figures
#   type    "joint" | "weighted" | "constrained" | "fg" | "csc" | "select"
#   link, baseline, df       for the direct models (directcif arguments)
#   cond    for "constrained": one G family, or several (then the one with the smallest AIC is used)
#   from    for "select": labels of earlier specifications; the one with the smallest AIC is used
#           (converged fits first). This is what an analyst would do without knowing the truth.
#   group   "direct (published)", "direct (SurvX)" or "reference", used for colours in the figures
model_formula <- Event(time, status) ~ x1 + x2

aic_of <- function(f) -2 * f@loglik + 2 * length(f@coefficients)
pick_best <- function(fits) {
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) return(NULL)
  ok <- vapply(fits, function(f) isTRUE(f@converged), logical(1))
  if (any(ok)) fits <- fits[ok]
  fits[[which.min(vapply(fits, aic_of, numeric(1)))]]
}
quietly <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)

fit_one <- function(spec, d) {
  switch(spec$type,
    joint = quietly(directcif(model_formula, d, baseline = if (is.null(spec$baseline)) "gompertz" else spec$baseline,
                              link = spec$link, df = if (is.null(spec$df)) 3 else spec$df)),
    weighted = quietly(directcif(model_formula, d, baseline = if (is.null(spec$baseline)) "gompertz" else spec$baseline,
                                 link = spec$link, estimation = "weighted", cause = 1)),
    constrained = pick_best(lapply(spec$cond, function(g)
      quietly(directcif(model_formula, d, link = spec$link, estimation = "constrained", cause = 1, cond = g)))),
    fg = quietly(fg(model_formula, d, cause = 1)),
    csc = quietly(csc(model_formula, d)),
    stop("unknown model type ", spec$type))
}

predict_one <- function(fit, spec, newdata, times) {
  if (is.null(fit)) return(NULL)
  P <- quietly(switch(spec$type,
    fg = predict(fit, newdata, times = times),
    csc = predict(fit, newdata, times = times, cause = 1),
    predict(fit, newdata, times = times, cause = 1)))
  if (is.null(P)) return(NULL)
  P <- matrix(as.numeric(P), nrow(newdata), length(times))
  if (any(!is.finite(P))) return(NULL)
  P
}

# fit every specification to the training data d; "select" entries refer to earlier fits
fit_models <- function(specs, d) {
  fits <- list(); secs <- numeric(0)
  for (s in specs) {
    t0 <- proc.time()[["elapsed"]]
    fits[s$label] <- list(if (s$type == "select") pick_best(fits[s$from]) else fit_one(s, d))
    secs[s$label] <- proc.time()[["elapsed"]] - t0
  }
  list(fits = fits, seconds = secs)
}

# -----------------------------------------------------------------------------------------------------
# 4. One replication, one cell, the whole study
# -----------------------------------------------------------------------------------------------------
# Every replication has its own seed, so results do not depend on the number of cores and a cell can
# be rerun alone. Data sets with too few events are redrawn (all models need
# events of both causes); the number of redraws is recorded.  (Rule: at least 5 events of cause 1 and
# 3 of cause 2.)
one_replication <- function(r, S, n, specs, test, seed) {
  set.seed(seed)
  redraws <- 0L
  repeat {
    d <- gen_data(n, S)
    if (sum(d$status == 1L) >= 5L && sum(d$status == 2L) >= 3L) break
    redraws <- redraws + 1L
    if (redraws > 100L) stop("scenario ", S$id, ", n = ", n, ": too few events")
  }
  fm <- fit_models(specs, d)
  rows <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    s <- specs[[i]]; f <- fm$fits[[s$label]]
    P <- predict_one(f, s, test$newdata, test$times)
    nt <- length(test$times)
    if (is.null(P)) {
      rows[[i]] <- data.frame(model = s$label, time = test$times, ok = FALSE, mse = NA_real_, bias = NA_real_, brier = NA_real_)
    } else {
      E <- P - test$truth
      rows[[i]] <- data.frame(model = s$label, time = test$times, ok = TRUE, mse = colMeans(E^2), bias = colMeans(E),
                              brier = colMeans((P - test$outcome)^2))
    }
    is_direct <- !is.null(f) && s$type %in% c("joint", "weighted", "constrained", "select")
    rows[[i]]$converged <- if (is.null(f)) NA else if (s$type %in% c("fg", "csc")) TRUE else isTRUE(f@converged)
    rows[[i]]$sum_above_1 <- if (is_direct && f@estimation == "joint") isTRUE(f@diagnostics$n_exceed_range > 0) else NA
    rows[[i]]$choice <- if (is_direct && s$type %in% c("select", "constrained"))
      paste0(f@estimation, ":", f@link[1L], if (f@estimation == "constrained") paste0(":", f@cond) else paste0(":", f@baseline[1L])) else NA_character_
    rows[[i]]$seconds <- fm$seconds[[s$label]]
  }
  out <- do.call(rbind, rows)
  cbind(scenario = S$id, n = n, rep = r, out,
        censored = mean(d$status == 0L), events1 = sum(d$status == 1L), events2 = sum(d$status == 2L), redraws = redraws)
}

# the test set: n_test independent subjects with their true risks and (uncensored) outcomes
make_test_set <- function(S, n_test, seed) {
  set.seed(seed)
  d <- gen_data(n_test, modifyList(S, list(cmax = Inf)))       # outcomes without censoring
  X <- as.matrix(d[, c("x1", "x2")])
  list(newdata = as.data.frame(X), times = S$times,
       truth = sapply(S$times, function(t) true_cif(X, t, S, 1L)),
       outcome = sapply(S$times, function(t) as.numeric(d$time <= t & d$status == 1L)))
}

# loads SurvX and this file on a parallel worker
source_on_worker <- function(f) { source(f); NULL }

# the function sent to the workers; built here so that it carries only what one replication needs
# (an error in one replication is recorded and does not stop the cell)
make_job <- function(S, n, specs, test, seeds) function(r)
  tryCatch(one_replication(r, S, n, specs, test, seeds[r]), error = function(e) conditionMessage(e))

run_cell <- function(S, n, R, specs, n_test, seed, cl = NULL) {
  si <- match(S$id, names(scenarios())); ni <- n
  test <- make_test_set(S, n_test, seed + 1000003L * si)
  seeds <- seed + 1000003L * si + 1009L * ni + seq_len(R)
  job <- make_job(S, n, specs, test, seeds)
  res <- if (is.null(cl)) lapply(seq_len(R), job) else parallel::parLapplyLB(cl, seq_len(R), job)
  failed <- !vapply(res, is.data.frame, logical(1))
  if (any(failed)) warning(sprintf("%s, n = %d: %d replication(s) failed, first error: %s", S$id, n, sum(failed),
                                   res[failed][[1]]), call. = FALSE)
  do.call(rbind, res[!failed])
}

# run all cells; every finished cell is saved to <out_dir>/cells/ and skipped when the study is
# restarted, so a long run can be interrupted and resumed. The combined raw results are written to
# <out_dir>/raw_<study>.csv.gz (one row per replication x model x horizon).
run_study <- function(study, specs, scenario_ids = names(scenarios()), sizes = c(100L, 250L, 500L, 1000L),
                      R = 1000L, n_test = 10000L, seed = 20261115L, cores = 1L, out_dir = "results", core_file) {
  dir.create(file.path(out_dir, "cells"), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(str(specs, max.level = 2)), file.path(out_dir, paste0("models_", study, ".txt")))
  cl <- NULL
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source_on_worker, core_file)
  }
  sc <- scenarios()
  t_all <- proc.time()[["elapsed"]]
  for (id in scenario_ids) for (n in sizes) {
    f <- file.path(out_dir, "cells", sprintf("%s_%s_n%04d_R%d.rds", study, id, n, R))
    if (file.exists(f)) { cat(sprintf("%-4s n = %4d  already done (%s)\n", id, n, basename(f))); next }
    t0 <- proc.time()[["elapsed"]]
    res <- run_cell(sc[[id]], n, R, specs, n_test, seed, cl)
    saveRDS(res, f)
    cat(sprintf("%-4s n = %4d  %4d replications  %6.0f s   RMSE x 100 (all horizons): %s\n", id, n, R,
                proc.time()[["elapsed"]] - t0,
                paste(sprintf("%s %.2f", unique(res$model),
                              tapply(res$mse, res$model, function(v) 100 * sqrt(mean(v, na.rm = TRUE)))[unique(res$model)]),
                      collapse = " | ")))
  }
  raw <- do.call(rbind, lapply(list.files(file.path(out_dir, "cells"), pattern = paste0("^", study, "_S.*_R", R, "\\.rds$"), full.names = TRUE), readRDS))
  con <- gzfile(file.path(out_dir, paste0("raw_", study, ".csv.gz")), "w"); utils::write.csv(raw, con, row.names = FALSE); close(con)
  cat(sprintf("\nall cells done in %.1f hours; raw results: %s\n", (proc.time()[["elapsed"]] - t_all) / 3600,
              file.path(out_dir, paste0("raw_", study, ".csv.gz"))))
  invisible(raw)
}

# -----------------------------------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------------------------------
# Per scenario x n x horizon x model:
#   rmse       sqrt of the mean over replications of the test-set MSE of the predicted risk of cause 1
#              against the true risk (Monte Carlo SE by the delta method)
#   rel_rmse   rmse / rmse of Fine-Gray, on the replications in which both models gave predictions,
#              with a 95% Monte Carlo interval (rel_lower, rel_upper);
#              reduction_pct = 100 (1 - rel_rmse) is the percentage reduction of the RMSE
#   win_pct    % of replications in which the model's MSE is below that of Fine-Gray
#   bias       mean signed error of the predicted risk, averaged over subjects and replications
#   dbrier     Brier score minus that of Fine-Gray (x 10^4), with its Monte Carlo z statistic
#   ok_pct, converged_pct, sum_above_1_pct, seconds
choice_text <- function(x) {
  if (all(is.na(x))) return(NA_character_)
  tb <- sort(table(x), decreasing = TRUE)
  paste(names(tb), as.integer(tb), collapse = "; ")
}

summarise_study <- function(raw, reference = "Fine-Gray") {
  key <- paste(raw$scenario, raw$n, raw$time, raw$rep)
  ref <- raw[raw$model == reference, ]
  rk <- paste(ref$scenario, ref$n, ref$time, ref$rep)
  raw$mse_ref <- ref$mse[match(key, rk)]; raw$brier_ref <- ref$brier[match(key, rk)]
  grp <- split(raw, list(raw$scenario, raw$n, raw$time, raw$model), drop = TRUE)
  out <- do.call(rbind, lapply(grp, function(g) {
    ok <- g$ok & is.finite(g$mse); pr <- ok & is.finite(g$mse_ref)
    m <- mean(g$mse[ok]); R <- sum(ok)
    a <- g$mse[pr]; b <- g$mse_ref[pr]; Rp <- length(a)
    rel <- sqrt(mean(a) / mean(b))
    # 95% Monte Carlo interval of rel by the delta method for the log of a ratio of paired means
    v <- if (Rp > 1L) stats::var(a) / (Rp * mean(a)^2) + stats::var(b) / (Rp * mean(b)^2) -
      2 * stats::cov(a, b) / (Rp * mean(a) * mean(b)) else NA_real_
    se_log <- 0.5 * sqrt(max(v, 0))
    db <- g$brier[pr] - g$brier_ref[pr]
    data.frame(scenario = g$scenario[1], n = g$n[1], time = g$time[1], model = g$model[1], replications = R,
               ok_pct = 100 * mean(g$ok), rmse = sqrt(m), rmse_mcse = stats::sd(g$mse[ok]) / sqrt(R) / (2 * sqrt(m)),
               rel_rmse = rel, rel_lower = rel * exp(-1.96 * se_log), rel_upper = rel * exp(1.96 * se_log),
               reduction_pct = 100 * (1 - rel),
               win_pct = if (g$model[1] == reference) NA else 100 * mean(g$mse[pr] < g$mse_ref[pr]),
               bias = mean(g$bias[ok]), bias_mcse = stats::sd(g$bias[ok]) / sqrt(R),
               brier = mean(g$brier[ok]), dbrier_x1e4 = if (g$model[1] == reference) NA else 1e4 * mean(db),
               dbrier_z = if (g$model[1] == reference || stats::sd(db) == 0) NA else mean(db) / (stats::sd(db) / sqrt(length(db))),
               converged_pct = 100 * mean(g$converged, na.rm = TRUE),
               sum_above_1_pct = if (all(is.na(g$sum_above_1))) NA else 100 * mean(g$sum_above_1, na.rm = TRUE),
               seconds = mean(g$seconds),
               choice = choice_text(g$choice))
  }))
  rownames(out) <- NULL
  out[order(out$scenario, out$n, out$time, out$model), ]
}

# design table: observed censoring and event counts per scenario and sample size
summarise_design <- function(raw) {
  one <- raw[!duplicated(paste(raw$scenario, raw$n, raw$rep)), ]
  agg <- stats::aggregate(cbind(censored_pct = 100 * censored, events1, events2, redraws) ~ scenario + n, one, mean)
  agg[order(agg$scenario, agg$n), ]
}

# locate the directory of the running script (Rscript file.R or source("file.R"))
this_dir <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  for (i in rev(seq_len(sys.nframe()))) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(dirname(normalizePath(o))) }
  getwd()
}
