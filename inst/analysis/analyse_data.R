#!/usr/bin/env Rscript
# Fit and compare the direct models, Fine-Gray and cause-specific Cox on one data set.
#
#   Rscript analyse_data.R                 runs on survival::mgus2 as a demonstration
#   edit the SETTINGS block for your own data, then run again
#
# Output (in the working directory): analysis_summaries.txt, analysis_coefficients.csv, analysis_fits.csv,
# analysis_comparison.csv, analysis_cif.pdf
suppressPackageStartupMessages(library(SurvX))

## ------------------------------- SETTINGS -------------------------------
data_file  <- NULL          # e.g. "mydata.csv"; NULL = demonstration data
time_var   <- "time"        # follow-up time
status_var <- "status"      # 0 = censored, 1, 2, ... = cause of failure
covariates <- NULL          # e.g. c("age", "sex", "stage"); NULL = all other columns
categorical <- NULL         # numeric-coded columns that are categories, e.g. c("stage", "centre")
cause      <- 1             # cause of interest
times      <- NULL          # prediction horizons; NULL = quartiles of the event times of the cause of interest
folds      <- 5             # cross-validation folds
seed       <- 1
## -------------------------------------------------------------------------

if (is.null(data_file)) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("set data_file, or install the survival package for the demonstration data")
  m <- survival::mgus2
  d <- data.frame(time = ifelse(m$pstat == 0, m$futime, m$ptime), status = ifelse(m$pstat == 0, 2L * m$death, 1L),
                  age = m$age, sex = m$sex, hgb = m$hgb, creat = m$creat, mspike = m$mspike)
  cat("Demonstration data: survival::mgus2; status 1 = progression to plasma cell malignancy, 2 = death without progression\n")
} else {
  d <- utils::read.csv(data_file, stringsAsFactors = TRUE)
}
if (is.null(covariates)) covariates <- setdiff(names(d), c(time_var, status_var))
stopifnot(all(c(time_var, status_var, covariates) %in% names(d)))
d <- d[, c(time_var, status_var, covariates)]
names(d)[1:2] <- c("time", "status")
for (v in categorical) d[[v]] <- factor(d[[v]])
for (v in covariates) if (is.character(d[[v]]) || is.logical(d[[v]])) d[[v]] <- factor(d[[v]])
n0 <- nrow(d); d <- d[stats::complete.cases(d) & d$time > 0, , drop = FALSE]
cat("rows:", n0, " used (complete, time > 0):", nrow(d), "\n")

## how each covariate enters the models
cat("\n== covariates ==\n")
for (v in covariates) {
  x <- d[[v]]
  if (is.factor(x)) cat(sprintf("  %-12s categorical, %d levels; reference = %s; indicator columns: %s\n", v, nlevels(x), levels(x)[1],
                                paste0(v, levels(x)[-1], collapse = ", ")))
  else cat(sprintf("  %-12s numeric, mean %.3g, sd %.3g (centred at its mean below)\n", v, mean(x), stats::sd(x)))
}
## numeric covariates are centred: for some baselines the baseline curve is the curve at covariate value 0
for (v in covariates) if (is.numeric(d[[v]])) d[[v]] <- d[[v]] - mean(d[[v]])

K <- max(d$status)
cat("\n== events ==\n"); print(table(status = d$status))
fml <- stats::as.formula(paste("Event(time, status) ~", paste(covariates, collapse = " + ")))
cat("\nmodel formula: "); print(fml)
if (is.null(times)) times <- as.numeric(stats::quantile(d$time[d$status == cause], c(0.25, 0.5, 0.75)))
cat("prediction horizons:", signif(times, 4), "\n")

## ------------------------------- models -------------------------------
quiet <- function(expr) tryCatch(expr, error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
spec <- list()
for (lk in c("PH", "PO", "GOR")) local({
  lk <- lk
  spec[[paste("direct joint", lk)]] <<- function(f, d) directcif(f, d, link = lk)
  spec[[paste("direct weighted", lk)]] <<- function(f, d) directcif(f, d, link = lk, estimation = "weighted", cause = cause)
})
if (K == 2) for (lk in c("PH", "PO")) for (cd in c("logistic", "weibull", "loglogistic", "lognormal")) local({
  lk <- lk; cd <- cd
  spec[[paste0("direct constrained ", lk, " (G ", cd, ")")]] <<- function(f, d) directcif(f, d, link = lk, estimation = "constrained", cause = cause, cond = cd)
})
spec[["Fine-Gray"]] <- function(f, d) fg(f, d, cause = cause)
spec[["cause-specific Cox"]] <- function(f, d) csc(f, d)

cat("\n== fitting", length(spec), "models to the full data ==\n")
fits <- lapply(spec, function(fun) quiet(fun(fml, d)))
fits <- fits[!vapply(fits, is.null, logical(1))]

## the full summary of every model, as the package prints it
summ <- unlist(lapply(names(fits), function(nm) c(paste0("################ ", nm, " ################"), capture.output(print(summary(fits[[nm]]))), "")))
writeLines(summ, "analysis_summaries.txt")
cat("summaries of all models: analysis_summaries.txt\n")
cat("\n== the main models ==\n")
for (nm in intersect(c("direct joint PH", "direct weighted PH", "Fine-Gray"), names(fits))) { cat("\n"); print(fits[[nm]]) }

is_direct <- function(f) S7::S7_inherits(f, SurvX:::DirectCIFFit)
fit_tab <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  if (is_direct(f)) data.frame(model = nm, converged = f@converged, loglik = f@loglik, npar = length(f@coefficients),
                               aic = if (f@estimation == "weighted") NA else AIC(f),
                               subjects_with_sum_of_CIFs_above_1 = f@diagnostics$n_exceed_range, notes = paste(SurvX:::.flag_text(f@flags), collapse = "; "))
  else data.frame(model = nm, converged = if (S7::S7_inherits(f, SurvX:::FineGrayFit)) f@converged else all(vapply(f@models, function(m) m$converged, logical(1))),
                  loglik = NA, npar = NA, aic = NA, subjects_with_sum_of_CIFs_above_1 = NA, notes = "")
}))
cat("\n== fit summary (AIC: comparable among the joint and constrained direct models only) ==\n")
print(fit_tab, row.names = FALSE, digits = 6)
utils::write.csv(fit_tab, "analysis_fits.csv", row.names = FALSE)

## coefficients of the cause of interest
coef_tab <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  if (is_direct(f)) {
    i <- f@index$cause == cause & f@index$role %in% c("beta", "alpha")
    if (!any(i)) return(NULL)
    data.frame(model = nm, term = f@index$term[i], estimate = unname(f@coefficients[i]), se = unname(sqrt(diag(f@vcov))[i]))
  } else if (S7::S7_inherits(f, SurvX:::FineGrayFit)) {
    data.frame(model = nm, term = names(f@coefficients), estimate = unname(f@coefficients), se = unname(sqrt(diag(f@vcov))))
  } else {
    m <- f@models[[cause]]
    data.frame(model = nm, term = names(m$coefficients), estimate = unname(m$coefficients), se = unname(sqrt(diag(m$vcov))))
  }
}))
coef_tab$z <- coef_tab$estimate / coef_tab$se; coef_tab$p <- 2 * stats::pnorm(-abs(coef_tab$z))
cat("\n== regression coefficients for cause", cause, "==\n",
    "direct PH and Fine-Gray: log subdistribution hazard ratios; direct PO: log odds ratios of the CIF;\n",
    "cause-specific Cox: log cause-specific hazard ratios (a different scale)\n", sep = "")
print(coef_tab, row.names = FALSE, digits = 4)
utils::write.csv(coef_tab, "analysis_coefficients.csv", row.names = FALSE)

## ------------------------------- cross-validated comparison (fitted and scored in the C++ core) -------------------------------
cat("\n== ", folds, "-fold cross-validation ==\n", sep = "")
direct <- c(lapply(c("PH", "PO", "GOR"), function(lk) c("gompertz", lk)), lapply(c("PH", "PO", "GOR"), function(lk) c("gompertz", lk, "weighted")))
if (K == 2) for (lk in c("PH", "PO")) for (cd in c("logistic", "weibull", "loglogistic", "lognormal")) direct <- c(direct, list(c("gompertz", lk, "constrained", cd)))
cmp <- compare_models(fml, d, times = times, direct = direct, cause = cause, folds = folds, seed = seed)
tab <- cmp@table
print(tab[order(tab$time, tab$brier_cv), c("time", "model", "brier_cv", "ipa_cv", "auc_cv", "calib_error_cv")], row.names = FALSE, digits = 4)
utils::write.csv(merge(tab, cmp@fit[, c("model", "aic", "folds_not_converged")]), "analysis_comparison.csv", row.names = FALSE)
cat("\nbest cross-validated Brier score at each horizon:\n")
for (t in unique(tab$time)) { s <- tab[tab$time == t & is.finite(tab$brier_cv), ]; cat(sprintf("  t = %-8.4g %s (Brier %.5f, IPA %.4f)\n", t, s$model[which.min(s$brier_cv)], min(s$brier_cv), s$ipa_cv[which.min(s$brier_cv)])) }

## ------------------------------- plot: average predicted CIF against Aalen-Johansen -------------------------------
grid <- seq(max(d$time) / 100, max(d$time), length.out = 100)
a <- aj(Event(time, status) ~ 1, d)@curves[[1]]
show <- c("direct joint PH", "direct weighted PH", grep("^direct constrained PH", names(fits), value = TRUE)[1], "Fine-Gray", "cause-specific Cox")
show <- intersect(show[!is.na(show)], names(fits))
grDevices::pdf("analysis_cif.pdf", width = 8, height = 6)
plot(c(0, a$time), c(0, a[[paste0("cif.", cause)]]), type = "s", lwd = 2, col = "grey40", xlab = "Time",
     ylab = paste("Cumulative incidence of cause", cause), main = "Average predicted CIF and the Aalen-Johansen estimate")
for (i in seq_along(show)) graphics::lines(grid, colMeans(predict(fits[[show[i]]], newdata = d, times = grid, cause = cause)), col = i + 1, lty = i)
graphics::legend("bottomright", legend = c("Aalen-Johansen", show), col = c("grey40", seq_along(show) + 1), lty = c(1, seq_along(show)), lwd = c(2, rep(1, length(show))), bty = "n")
invisible(grDevices::dev.off())
cat("\nwritten: analysis_summaries.txt, analysis_fits.csv, analysis_coefficients.csv, analysis_comparison.csv, analysis_cif.pdf\n")
