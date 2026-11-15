#!/usr/bin/env Rscript
# Checks the R layer of SurvX against survival, cmprsk and riskRegression.
#   Rscript check_SurvX.R        (after R CMD INSTALL SurvX)
# Writes check_SurvX_report.csv in the working directory -- please send it back.
# Written without access to R. The checks up to version 0.3.1 were run by the user (all passed); the
# block "output methods" is new (versions 0.4.x) and has not been run. Every block is wrapped so that an error is recorded and the script continues.
suppressPackageStartupMessages({ library(SurvX); library(survival) })
if (packageVersion("SurvX") < "0.5.1")
  stop("this script needs SurvX 0.5.1 or later, but version ", packageVersion("SurvX"), " is loaded from ",
       find.package("SurvX"), ".\nInstall the new tar.gz in a fresh R session (restart R first), then run the script again.")
cat("SurvX", format(packageVersion("SurvX")), "from", find.package("SurvX"), "\n")
has_cmprsk <- requireNamespace("cmprsk", quietly = TRUE)
has_rr <- requireNamespace("riskRegression", quietly = TRUE) && requireNamespace("prodlim", quietly = TRUE)
if (has_rr) Hist <- prodlim::Hist     # riskRegression::Score recognises the response by the name Hist
report <- NULL
add <- function(item, reference, ours, ref, tol = 1e-6) {
  ours <- as.numeric(ours); ref <- as.numeric(ref)
  if (length(ours) != length(ref) || !length(ours)) { st <- sprintf("LENGTH %d vs %d", length(ours), length(ref)); ad <- NA }
  else { ok <- is.finite(ours) & is.finite(ref); ad <- if (any(ok)) max(abs(ours[ok] - ref[ok]) / pmax(1, abs(ref[ok]))) else NA
         st <- if (any(is.finite(ours) != is.finite(ref))) "DIFF (non-finite)" else if (isTRUE(ad <= tol)) "PASS" else "DIFF" }
  report <<- rbind(report, data.frame(item, reference, max_scaled_diff = ad, tol, status = st))
}
block <- function(item, expr) tryCatch(expr, error = function(e)
  report <<- rbind(report, data.frame(item, reference = "", max_scaled_diff = NA, tol = NA, status = paste("ERROR:", conditionMessage(e)))))

m <- survival::mgus2
m$time <- ifelse(m$pstat == 0, m$futime, m$ptime); m$status <- ifelse(m$pstat == 0, 2L * m$death, 1L)
d <- na.omit(data.frame(time = m$time, status = m$status, dead = as.integer(m$status > 0), sex = m$sex,
                        age = m$age - 70, hgb = m$hgb - 13, mspike = m$mspike))
d <- d[d$time > 0, ]

block("Event / print methods", { y <- Event(d$time, d$status); stopifnot(inherits(y[1:5], "Event")); capture.output(print(y[1:5])); add("Event subsetting", "base R", nrow(unclass(y[1:5])), 5) })

block("km", {
  k <- km(Event(time, dead) ~ sex, d); capture.output(print(k)); capture.output(summary(k, times = c(50, 100)))
  s <- survfit(Surv(time, dead) ~ sex, d); ss <- summary(s, censored = FALSE)
  ours <- do.call(rbind, k@curves); ours <- ours[ours$n.event > 0, ]
  add("km: survival", "survival::survfit", ours$surv, ss$surv)
  add("km: Greenwood se", "survival::survfit", ours$se, ss$std.err)
  sd <- survdiff(Surv(time, dead) ~ sex, d)
  add("km: log-rank chi-square", "survival::survdiff", k@test$statistic, sd$chisq)
  pdf(NULL); plot(k); dev.off()
})
block("rmst", { s <- survfit(Surv(time, dead) ~ 1, d); tb <- summary(s, rmean = 200)$table
  add("rmst", "survival (rmean)", rmst(d$time, d$dead, 200), c(tb[["rmean"]], tb[["se(rmean)"]]), 1e-6) })

block("aj", {
  a <- aj(Event(time, status) ~ sex, d); capture.output(print(a)); pdf(NULL); plot(a, cause = 1); dev.off()
  if (has_cmprsk) {
    ci <- cmprsk::cuminc(d$time, d$status, d$sex)
    tp <- c(50, 100, 200); tt <- cmprsk::timepoints(ci, tp)
    ours <- summary(a, times = tp)@tables
    add("aj: CIF of cause 1 at 50/100/200, group F", "cmprsk::cuminc", ours[["F"]]$cif.1, tt$est["F 1", ])
    add("aj: CIF of cause 2 at 50/100/200, group M", "cmprsk::cuminc", ours[["M"]]$cif.2, tt$est["M 2", ])
    add("aj: Gray test statistics", "cmprsk::cuminc", a@tests$statistic, ci$Tests[, "stat"])
  }
})

block("cox", {
  f <- cox(Event(time, dead) ~ sex + age + hgb + mspike, d); capture.output(print(f)); r <- coxph(Surv(time, dead) ~ sex + age + hgb + mspike, d)
  add("cox: coefficients", "survival::coxph", coef(f), coef(r)); add("cox: standard errors", "survival::coxph", sqrt(diag(vcov(f))), sqrt(diag(vcov(r))))
  add("cox: log-likelihood", "survival::coxph", as.numeric(logLik(f)), r$loglik[2])
  nd <- d[1:5, ]; s <- summary(survfit(r, newdata = nd), times = c(50, 100))$surv
  add("cox: predicted survival", "survival::survfit.coxph", predict(f, nd, times = c(50, 100), type = "survival"), t(s), 1e-6)
  fb <- cox(Event(time, dead) ~ age, d, ties = "breslow"); rb <- coxph(Surv(time, dead) ~ age, d, ties = "breslow")
  add("cox (Breslow): coefficient", "survival::coxph", coef(fb), coef(rb))
})

block("fg", {
  f <- fg(Event(time, status) ~ sex + age + hgb + mspike, d); capture.output(print(f))
  if (has_cmprsk) {
    X <- model.matrix(~ sex + age + hgb + mspike, d)[, -1]
    r <- cmprsk::crr(d$time, d$status, X, failcode = 1, cencode = 0, gtol = 1e-10)
    add("fg: coefficients", "cmprsk::crr", coef(f), r$coef); add("fg: robust standard errors", "cmprsk::crr", sqrt(diag(vcov(f))), sqrt(diag(r$var)), 1e-5)
    pr <- predict(r, X[1:5, ]); tp <- c(50, 100)
    refp <- sapply(tp, function(s) pr[max(which(pr[, 1] <= s)), -1])
    add("fg: predicted CIF", "cmprsk::predict.crr", predict(f, d[1:5, ], times = tp), refp, 1e-6)
  }
})

block("csc", {
  f <- csc(Event(time, status) ~ sex + age + hgb + mspike, d); capture.output(print(f))
  r1 <- coxph(Surv(time, status == 1) ~ sex + age + hgb + mspike, d)
  add("csc: coefficients of cause 1", "survival::coxph", coef(f)[[1]], coef(r1))
  if (has_rr) {
    rr <- riskRegression::CSC(Hist(time, status) ~ sex + age + hgb + mspike, data = d)
    add("csc: predicted CIF of cause 1", "riskRegression::predictRisk", predict(f, d[1:5, ], times = c(50, 100), cause = 1),
        riskRegression::predictRisk(rr, newdata = d[1:5, ], times = c(50, 100), cause = 1), 1e-6)
  }
})

# independent R code of the direct Gompertz likelihood (Jeong & Fine 2007, eq. 3.9, 4.1-4.2)
refll <- function(th, time, status, X, alpha) {
  p <- ncol(X); S <- rep(1, length(time)); ll <- 0
  for (k in 1:2) {
    q <- th[(k - 1) * (2 + p) + 1:(2 + p)]; rho <- q[1]; tau <- exp(q[2]); eta <- drop(X %*% q[-(1:2)])
    u <- tau * expm1(rho * time) / rho; du <- tau * exp(rho * time); eu <- exp(eta) * u
    if (alpha == 0) { Fk <- -expm1(-eu); lf <- eta + log(du) - eu } else { Fk <- 1 - (1 + alpha * eu)^(-1 / alpha); lf <- eta + log(du) - (1 / alpha + 1) * log1p(alpha * eu) }
    ll <- ll + sum(lf[status == k]); S <- S - Fk
  }
  ll + sum(log(S[status == 0]))
}
block("directcif", {
  X <- model.matrix(~ sex + age + hgb + mspike, d)[, -1]
  for (lk in c("PH", "PO")) {
    f <- directcif(Event(time, status) ~ sex + age + hgb + mspike, d, baseline = "gompertz", link = lk); capture.output(print(f))
    add(paste("directcif Gompertz", lk, ": log-likelihood"), "R code of the paper formulas", as.numeric(logLik(f)), refll(coef(f), d$time, d$status, X, if (lk == "PH") 0 else 1), 1e-9)
    add(paste("directcif Gompertz", lk, ": AIC"), "-2 logLik + 2 npar", AIC(f), -2 * as.numeric(logLik(f)) + 2 * length(coef(f)), 1e-12)
    pr <- predict(f, d[1:5, ], times = c(50, 100, Inf), cause = 1, se = TRUE)
    q <- coef(f)[1:6]; eu <- exp(drop(X[1:5, ] %*% q[3:6])) %o% (exp(q[2]) * ifelse(is.finite(c(50, 100, Inf)), expm1(q[1] * c(50, 100, Inf)) / q[1], -1 / q[1]))
    add(paste("directcif Gompertz", lk, ": predicted CIF and plateau"), "R code of eq. (3.9), (3.10)", pr$cif, if (lk == "PH") -expm1(-eu) else 1 - 1 / (1 + eu), 1e-9)
    add(paste("directcif Gompertz", lk, ": prediction se finite and positive"), "delta method", as.numeric(all(is.finite(pr$se) & pr$se > 0)), 1)
  }
  for (b in c("gompertz3", "weibull", "mix_weibull", "mix_loglogistic", "mix_lognormal", "logistic3", "spline")) {
    f <- directcif(Event(time, status) ~ sex + age + hgb + mspike, d, baseline = b); capture.output(print(f))
    pr <- predict(f, d[1:5, ], times = c(50, 100))
    add(paste("directcif", b, ": runs, predictions in [0, 1]"), "sanity", as.numeric(all(pr >= 0 & pr <= 1)), 1)
  }
  f0 <- directcif(Event(time, status) ~ 1, d); capture.output(print(f0))
  add("directcif without covariates (Jeong & Fine 2006): 4 parameters", "sanity", length(coef(f0)), 4)
  fg0 <- directcif(Event(time, status) ~ age, d, link = "GOR"); capture.output(print(fg0))
  pdf(NULL); plot(f0); plot(directcif(Event(time, status) ~ sex + age, d), newdata = data.frame(sex = c("F", "M"), age = 0)); dev.off()
})

block("metrics and comparison", {
  f <- directcif(Event(time, status) ~ sex + age + hgb + mspike, d); risk <- drop(predict(f, times = 100))
  b <- brier_score(d$time, d$status, risk, 100); a <- auc_t(d$time, d$status, risk, 100); cal <- calibration(d$time, d$status, risk, 100)
  if (has_rr) {
    sc <- riskRegression::Score(list(m = matrix(risk)), formula = Hist(time, status) ~ 1, data = d, times = 100, cause = 1,
                                metrics = c("auc", "brier"), summary = "ipa", cens.model = "km", conf.int = FALSE)
    bs <- as.data.frame(sc$Brier$score); au <- as.data.frame(sc$AUC$score)
    add("brier_score", "riskRegression::Score", b[["brier"]], bs$Brier[bs$model == "m"], 1e-6)
    add("auc_t", "riskRegression::Score", a, au$AUC[au$model == "m"], 1e-6)
  }
  add("calibration: groups", "sanity", nrow(cal), 10)
  cm <- compare_models(Event(time, status) ~ sex + age + hgb + mspike, d, times = c(50, 100), folds = 5)
  print(cm)
  add("compare_models: rows (default list: 5 direct models, Fine-Gray, cause-specific Cox)", "sanity", nrow(cm@table), 2 * 7)
  fwc <- directcif(Event(time, status) ~ sex + age + hgb + mspike, d, estimation = "weighted")
  add("compare_models: apparent Brier of the weighted direct PH row at t = 100", "brier_score()", cm@table$brier[cm@table$model == "direct gompertz PH, weighted" & cm@table$time == 100],
      brier_score(d$time, d$status, drop(predict(fwc, times = 100)), 100)[["brier"]], 1e-9)
  fcc <- directcif(Event(time, status) ~ sex + age + hgb + mspike, d, estimation = "constrained", cond = "weibull")
  add("compare_models: apparent Brier of the constrained direct PH row at t = 100", "brier_score()", cm@table$brier[cm@table$model == "direct gompertz PH, constrained (G weibull)" & cm@table$time == 100],
      brier_score(d$time, d$status, drop(predict(fcc, times = 100)), 100)[["brier"]], 1e-9)
  add("compare_models: AIC reported for joint and constrained rows only", "sanity", as.numeric(is.na(cm@fit$aic)), c(0, 0, 1, 1, 0, 1, 1))
  cm2 <- compare_models(Event(time, status) ~ sex + age, d, times = 100, folds = 3, direct = list(c("weibull", "PO", "weighted"), c("logistic3", "PH", "constrained", "lognormal"), c("spline", "PH")))
  add("compare_models with a user list: rows", "sanity", nrow(cm2@table), 5)
  add("compare_models: apparent Brier of the direct Gompertz PH row at t = 100", "brier_score()", cm@table$brier[cm@table$model == "direct gompertz PH" & cm@table$time == 100], b[["brier"]], 1e-9)
})

# ---- estimation = "weighted" and "constrained": independent R code of the likelihoods ----
ref_weighted <- function(th, time, status, X, alpha, cause = 1) {
  rho <- th[1]; tau <- exp(th[2]); eta <- drop(X %*% th[-(1:2)])
  eu <- function(t, e) exp(e) * tau * expm1(rho * t) / rho
  lS <- function(t, e) if (alpha == 0) -eu(t, e) else -log1p(alpha * eu(t, e)) / alpha
  lf <- function(t, e) e + log(tau) + rho * t - (if (alpha == 0) eu(t, e) else (1 / alpha + 1) * log1p(alpha * eu(t, e)))
  ct <- sort(unique(time[status == 0])); G <- 1; mass <- numeric(length(ct))        # Kaplan-Meier of the censoring times
  for (j in seq_along(ct)) { dc <- sum(time == ct[j] & status == 0); r <- sum(time > ct[j]) + dc; Gn <- G * (1 - dc / r); mass[j] <- G - Gn; G <- Gn }
  if (G > 1e-14) { ct <- c(ct, max(time)); mass <- c(mass, G) }
  ll <- 0
  for (i in seq_along(time)) {
    if (status[i] == cause) ll <- ll + lf(time[i], eta[i])
    else if (status[i] == 0) ll <- ll + lS(time[i], eta[i])
    else { k <- ct >= time[i]; if (any(k)) ll <- ll + sum(mass[k] / sum(mass[k]) * lS(ct[k], eta[i])) }
  }
  ll
}
ref_constrained <- function(th, time, status, X) {       # Gompertz PH for cause 1, Weibull G (Shi et al. 2013, eq. 5, 7)
  p <- ncol(X); rho <- th[1]; tau <- exp(th[2]); eta <- drop(X %*% th[2 + seq_len(p)]); lam <- exp(th[p + 3]); ka <- exp(th[p + 4])
  Ft <- -expm1(-exp(eta) * tau * expm1(rho * time) / rho); Finf <- -expm1(exp(eta) * tau / rho)
  lf1 <- eta + log(tau) + rho * time - exp(eta) * tau * expm1(rho * time) / rho
  G <- pweibull(time, ka, 1 / lam); g <- dweibull(time, ka, 1 / lam)
  sum(lf1[status == 1]) + sum(log(1 - Finf[status == 2]) + log(g[status == 2])) + sum(log(1 - Ft - (1 - Finf) * G)[status == 0])
}
block("directcif: weighted and constrained", {
  X <- model.matrix(~ sex + age + hgb + mspike, d)[, -1]; fm <- Event(time, status) ~ sex + age + hgb + mspike
  for (lk in c("PH", "PO")) {
    fw <- directcif(fm, d, link = lk, estimation = "weighted", cause = 1); capture.output(print(fw))
    add(paste("weighted", lk, ": log pseudo-likelihood"), "R code of the weighted likelihood", as.numeric(logLik(fw)), ref_weighted(coef(fw), d$time, d$status, X, if (lk == "PH") 0 else 1), 1e-9)
    add(paste("weighted", lk, ": converged"), "sanity", as.numeric(fw@converged), 1)
  }
  fw <- directcif(fm, d, estimation = "weighted", cause = 1); ff <- fg(fm, d)
  add("weighted PH: regression coefficients close to Fine-Gray (same estimand, different baseline)", "SurvX::fg", coef(fw)[-(1:2)], coef(ff), 0.1)
  add("weighted PH: sandwich se close to Fine-Gray robust se", "SurvX::fg", sqrt(diag(vcov(fw)))[-(1:2)], sqrt(diag(vcov(ff))), 0.1)
  pw <- predict(fw, d[1:5, ], times = c(50, 100, Inf), se = TRUE)
  add("weighted: predictions in [0, 1], se finite", "sanity", as.numeric(all(pw$cif >= 0 & pw$cif <= 1 & is.finite(pw$se))), 1)
  f2 <- directcif(fm, d, estimation = "weighted", cause = 2, baseline = "weibull"); fs <- directcif(fm, d, estimation = "weighted", baseline = "spline")
  add("weighted for cause 2 (Weibull) and with a spline baseline: run", "sanity", as.numeric(is.finite(logLik(f2)) && is.finite(logLik(fs))), 1)
  fc <- directcif(fm, d, estimation = "constrained", cause = 1, cond = "weibull"); capture.output(print(fc))
  add("constrained (Gompertz PH, Weibull G): log-likelihood", "R code of Shi et al. (2013), eq. (5), (7)", as.numeric(logLik(fc)), ref_constrained(coef(fc), d$time, d$status, X), 1e-9)
  p1 <- predict(fc, d[1:5, ], times = c(50, Inf), cause = 1); p2 <- predict(fc, d[1:5, ], times = c(50, Inf), cause = 2, se = TRUE)
  add("constrained: asymptotes of the two causes add up to one", "Shi et al. (2013), eq. (1)", p1[, 2] + p2$cif[, 2], rep(1, 5), 1e-12)
  add("constrained: sum of the CIFs below one at t = 50, se finite", "sanity", as.numeric(all(p1[, 1] + p2$cif[, 1] < 1) && all(is.finite(p2$se))), 1)
  for (cd in c("logistic", "loglogistic", "lognormal")) {
    f <- directcif(fm, d, link = "PO", estimation = "constrained", cause = 1, cond = cd)
    add(paste("constrained PO, G", cd, ": runs"), "sanity", as.numeric(is.finite(logLik(f))), 1)
  }
  cc <- cv_compare(fm, d, times = c(50, 100), folds = 3, models = list(
    "weighted PH" = function(f, d) directcif(f, d, estimation = "weighted"),
    "constrained PH" = function(f, d) directcif(f, d, estimation = "constrained", cond = "weibull"),
    "Fine-Gray" = function(f, d) fg(f, d), "cause-specific Cox" = function(f, d) csc(f, d)))
  print(cc)
  add("cv_compare: rows", "sanity", nrow(cc@table), 8)
  add("cv_compare: apparent Brier of Fine-Gray at t = 100", "brier_score()", cc@table$brier[cc@table$model == "Fine-Gray" & cc@table$time == 100],
      brier_score(d$time, d$status, drop(predict(ff, times = 100)), 100)[["brier"]], 1e-9)
  add("cv_compare: cross-validated measures finite", "sanity", as.numeric(all(is.finite(cc@table$brier_cv))), 1)
})

# ---- factors and character covariates ----
block("categorical covariates", {
  d2 <- d; d2$grp <- c("a", "b", "c")[1 + (seq_len(nrow(d2)) %% 3)]; d2$sexn <- as.integer(d2$sex == "M")
  f1 <- fg(Event(time, status) ~ grp + age, d2); f2 <- fg(Event(time, status) ~ factor(sexn) + age, d2); f3 <- fg(Event(time, status) ~ sex + age, d2)
  add("character column becomes indicator columns (2 for 3 levels)", "sanity", length(coef(f1)), 3)
  add("factor(numeric code) = factor column", "SurvX::fg", unname(coef(f2)), unname(coef(f3)), 1e-12)
  pr <- predict(f1, data.frame(grp = c("c", "a"), age = 0), times = 100)
  add("prediction for new data given as level names", "sanity", as.numeric(all(is.finite(pr)) && nrow(pr) == 2), 1)
})

# ---- user-facing output: every print and summary goes to SurvX_output_examples.txt, every plot to SurvX_plots.pdf ----
block("output methods", {
  fm <- Event(time, status) ~ sex + age + hgb + mspike
  k <- km(Event(time, dead) ~ sex, d); a <- aj(Event(time, status) ~ sex, d)
  fc <- cox(Event(time, dead) ~ sex + age + hgb + mspike, d); ff <- fg(fm, d); fs <- csc(fm, d)
  fj <- directcif(fm, d); fpo <- update(fj, link = "PO"); fgo <- update(fj, link = "GOR")
  fw <- directcif(fm, d, estimation = "weighted"); fk <- directcif(fm, d, estimation = "constrained", cond = "weibull")
  fsp <- directcif(fm, d, baseline = "spline")
  cm <- compare_models(fm, d, times = c(50, 100), folds = 3)
  txt <- c(paste("SurvX", format(packageVersion("SurvX")), "| script for 0.5.1 | R", getRversion()), "",
    "######## km ########", capture.output(print(k)), "", capture.output(print(summary(k, times = c(50, 100, 200)))),
    "", "######## aj ########", capture.output(print(a)), "", capture.output(print(summary(a, times = c(50, 100, 200)))),
    "", "######## cox ########", capture.output(print(fc)), "", capture.output(print(summary(fc))),
    "", "######## fg ########", capture.output(print(ff)), "", capture.output(print(summary(ff))),
    "", "######## csc ########", capture.output(print(fs)), "", capture.output(print(summary(fs))),
    "", "######## directcif joint PH ########", capture.output(print(fj)), "", capture.output(print(summary(fj))),
    "", "######## directcif joint PO ########", capture.output(print(fpo)),
    "", "######## directcif joint GOR ########", capture.output(print(fgo)),
    "", "######## directcif spline ########", capture.output(print(fsp)),
    "", "######## directcif weighted ########", capture.output(print(fw)), "", capture.output(print(summary(fw))),
    "", "######## directcif constrained ########", capture.output(print(fk)), "", capture.output(print(summary(fk))),
    "", "######## estimates(fj) ########", capture.output(print(estimates(fj))),
    "", "######## confint(fj) ########", capture.output(print(confint(fj))),
    "", "######## AIC(fj, fpo, fgo) ########", capture.output(print(AIC(fj, fpo, fgo))),
    "", "######## compare_models ########", capture.output(print(cm)))
  writeLines(txt, "SurvX_output_examples.txt")
  grDevices::pdf("SurvX_plots.pdf", width = 10, height = 6.5)
  graphics::plot.new(); graphics::text(0.5, 0.5, paste("SurvX", format(packageVersion("SurvX")), "- plots of check_SurvX.R for 0.5.1"), cex = 1.5)
  nd <- data.frame(sex = c("F", "M"), age = 0, hgb = 0, mspike = 1)
  plot(k, main = "km(): two groups"); plot(km(Event(time, dead) ~ 1, d), main = "km(): one group, median", median = TRUE)
  plot(k, conf.int = TRUE, risk.table = FALSE, percent = TRUE, mark.times = c(100, 200), main = "km(): options")
  plot(a, main = "aj(): two groups"); plot(aj(Event(time, status) ~ 1, d), layout = "stacked", main = "aj(): stacked")
  plot(fc, newdata = nd, labels = c("female", "male"))
  cif_plot(fj, main = "direct joint PH: typical covariate values")
  cif_plot(fj, newdata = nd, labels = c("female", "male"), main = "direct joint PH: two profiles")
  cif_plot(fj, type = "marginal", reference = TRUE, main = "direct joint PH: average curve and Aalen-Johansen")
  cif_plot(fj, type = "unadjusted", reference = TRUE, main = "direct joint PH without covariates and Aalen-Johansen")
  cif_plot(fj, layout = "overlay", main = "overlay of the causes"); cif_plot(fj, newdata = nd, layout = "overlay", main = "overlay, two profiles")
  cif_plot(fj, newdata = nd, layout = "stacked", main = "stacked, two profiles")
  cif_plot(list("direct PH" = fj, "direct PO" = fpo, "weighted" = fw, "constrained" = fk, "Fine-Gray" = ff, "cause-specific Cox" = fs),
           cause = 1, B = 100, main = "six models, cause 1")
  cif_plot(list("direct PH" = fj, "Fine-Gray" = ff, "cause-specific Cox" = fs), B = 100, cause.names = c("progression", "death"),
           percent = TRUE, grid = TRUE, mark.times = c(120, 240), main = "three models, all causes, options")
  cif_plot(list("direct PH" = fj, "Fine-Gray" = ff), conf.int = FALSE, risk.table = FALSE, legend = FALSE, caption = FALSE, main = "everything switched off")
  plot(fw, newdata = nd); plot(fk, newdata = nd, cause = 2); plot(fs, B = 100); plot(ff, B = 100)
  plot(cm); plot(cm, measure = "brier"); plot(cm, measure = "auc", highlight = FALSE, main = "AUC of all models")
  cif_plot(fj, main = paste("a deliberately long title to check that titles and captions are wrapped",
                            "to the width of the page instead of running off it, also on a narrow device or with a large cex"))
  plot(k, main = paste("Kaplan-Meier curves by sex with a deliberately long title that has to be wrapped",
                       "onto a second line so that it stays inside the page"))
  invisible(grDevices::dev.off())
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    grDevices::pdf("SurvX_plots_ggplot2.pdf", width = 10, height = 6.5)
    print(ggplot2::autoplot(fj)); print(cif_plot(list("direct PH" = fj, "Fine-Gray" = ff), B = 50, engine = "ggplot2"))
    print(cif_plot(fj, newdata = nd, layout = "stacked", engine = "ggplot2"))
    invisible(grDevices::dev.off())
  }
  # cif_data() against predict()
  mode <- names(which.max(table(d$sex)))
  prof <- data.frame(sex = factor(mode, levels = levels(d$sex)), age = mean(d$age), hgb = mean(d$hgb), mspike = mean(d$mspike))
  tp <- c(50, 100, 200)
  cd <- function(...) { r <- cif_data(..., times = tp); r }
  add("cif_data: typical values = mean of numeric, most frequent level (direct, both causes)", "predict()",
      cd(fj, conf.int = FALSE)$cif, c(predict(fj, prof, times = tp, cause = 1), predict(fj, prof, times = tp, cause = 2)), 1e-10)
  add("cif_data: Fine-Gray at typical values, cause 1", "predict()", cd(ff, cause = 1, conf.int = FALSE)$cif, predict(ff, prof, times = tp), 1e-10)
  add("cif_data: Fine-Gray refitted for cause 2", "fg(cause = 2)", cd(ff, cause = 2, conf.int = FALSE)$cif, predict(fg(fm, d, cause = 2), prof, times = tp), 1e-8)
  add("cif_data: cause-specific Cox at typical values, both causes", "predict()", cd(fs, conf.int = FALSE)$cif,
      c(predict(fs, prof, times = tp, cause = 1), predict(fs, prof, times = tp, cause = 2)), 1e-10)
  add("cif_data: weighted model refitted for cause 2", "directcif(estimation = 'weighted', cause = 2)", cd(fw, cause = 2, conf.int = FALSE)$cif,
      predict(directcif(fm, d, estimation = "weighted", cause = 2), prof, times = tp), 1e-8)
  add("cif_data: constrained model, cause 2", "predict()", cd(fk, cause = 2, conf.int = FALSE)$cif, predict(fk, prof, times = tp, cause = 2), 1e-10)
  add("cif_data: marginal = mean of the predicted curves (direct)", "predict()", cd(fj, type = "marginal", conf.int = FALSE, cause = 1)$cif,
      colMeans(predict(fj, times = tp, cause = 1)), 1e-10)
  add("cif_data: marginal = mean of the predicted curves (Fine-Gray)", "predict()", cd(ff, type = "marginal", conf.int = FALSE, cause = 1)$cif,
      colMeans(predict(ff, times = tp)), 1e-10)
  add("cif_data: direct model without covariates", "directcif(Event(time, status) ~ 1)", cd(fj, type = "unadjusted", conf.int = FALSE, cause = 1)$cif,
      drop(predict(directcif(Event(time, status) ~ 1, d), data.frame(row = 1), times = tp, cause = 1)), 1e-8)
  ajr <- summary(aj(Event(time, status) ~ 1, d), times = tp)@tables[[1]]
  add("cif_data: cause-specific Cox without covariates = Aalen-Johansen", "aj()", cd(fs, type = "unadjusted", conf.int = FALSE)$cif, c(ajr$cif.1, ajr$cif.2), 1e-10)
  add("cif_data: Fine-Gray without covariates close to Aalen-Johansen", "aj()", cd(ff, type = "unadjusted", conf.int = FALSE, cause = 1)$cif, ajr$cif.1, 0.01)
  b1 <- cd(ff, cause = 1, B = 50, seed = 7); b2 <- cd(ff, cause = 1, B = 50, seed = 7)
  add("cif_data: bootstrap reproducible with the same seed", "sanity", c(b1$lower, b1$upper), c(b2$lower, b2$upper), 0)
  add("cif_data: bootstrap limits around the estimate", "sanity", as.numeric(all(b1$lower <= b1$cif & b1$cif <= b1$upper)), 1)
  dj <- cd(fj, newdata = prof, cause = 1); pj <- predict(fj, prof, times = tp, cause = 1, se = TRUE)
  add("cif_data: delta-method se of the direct model", "predict(se = TRUE)", dj$se, drop(pj$se), 1e-10)

  ss <- summary(survfit(Surv(time, dead) ~ sex, d), censored = FALSE); ours <- do.call(rbind, k@curves); ours <- ours[ours$n.event > 0, ]
  add("km: lower and upper 95% limits (log scale)", "survival::survfit", c(ours$lower, ours$upper), c(ss$lower, ss$upper), 1e-6)
  med <- summary(survfit(Surv(time, dead) ~ sex, d))$table
  add("km: median survival and its limits", "survival::survfit", c(k@counts$median, k@counts$median.lower, k@counts$median.upper),
      c(med[, "median"], med[, "0.95LCL"], med[, "0.95UCL"]), 1e-9)
  e <- estimates(fj); ci <- confint(fj)
  add("estimates() and confint() agree", "confint", c(e$lower, e$upper), c(ci[, 1], ci[, 2]), 1e-12)
  ef <- estimates(ff)
  add("estimates(fg): exponentiated effects", "exp(coef)", ef$exp_estimate, exp(coef(ff)), 1e-12)
  add("update(): link changed, data kept", "sanity", as.numeric(c(fpo@link[1] == "PO", nobs(fpo) == nobs(fj))), c(1, 1))
  add("AIC() of three direct fits", "-2 logLik + 2 npar", AIC(fj, fpo, fgo)$AIC, sapply(list(fj, fpo, fgo), function(f) -2 * as.numeric(logLik(f)) + 2 * length(coef(f))), 1e-9)
  add("summary() keeps the coefficient table", "estimates()", as.data.frame(summary(fj))$estimate, estimates(fj)$estimate, 0)
  add("model.matrix() and formula()", "sanity", as.numeric(c(ncol(model.matrix(fj)) == 4, inherits(formula(fj), "formula"))), c(1, 1))
})

options(width = 200)
print(report, row.names = FALSE, digits = 3)
cat("\nstatus counts:\n"); print(table(sub(":.*", "", report$status)))
write.csv(report, "check_SurvX_report.csv", row.names = FALSE)
cat("\nR", as.character(getRversion()), "| S7", as.character(packageVersion("S7")), "| SurvX", as.character(packageVersion("SurvX")),
    "\nwritten: check_SurvX_report.csv, SurvX_output_examples.txt, SurvX_plots.pdf (and SurvX_plots_ggplot2.pdf with ggplot2)\n")
