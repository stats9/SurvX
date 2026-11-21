# Paper 2: real-data analyses. Run from the repository root:
#   Rscript paper2/realdata/realdata.R
# Writes paper2/results/realdata.txt, CSV tables and PDF figures.
library(SurvX)
out_dir <- file.path("paper2", "results")
dir.create(out_dir, showWarnings = FALSE)
safe <- function(expr) tryCatch(expr, error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
se_tab <- function(fit) {
  e <- estimates(fit)
  data.frame(term = paste0(e$cause, ":", e$term), estimate = e$estimate,
             se_ipcw = sqrt(diag(fit@diagnostics$vcov_ipcw)), se_ignored = sqrt(diag(fit@diagnostics$vcov_robust)),
             ratio = sqrt(diag(fit@diagnostics$vcov_ipcw)) / sqrt(diag(fit@diagnostics$vcov_robust)))
}

sink(file.path(out_dir, "realdata.txt"), split = TRUE)
cat("SurvX", as.character(utils::packageVersion("SurvX")), "\n")

# ------------------------------------------------------------------------------------------
# 1. fertility: two causes (live birth; stillbirth or abortion). Weighted estimator for each
#    cause, links, and the standard errors with and without the censoring term.
# ------------------------------------------------------------------------------------------
cat("\n==== fertility ====\n")
data(fertility, package = "SurvX")
fertility$age_c <- fertility$age - round(mean(fertility$age))
print(table(fertility$status))
print(aj(Event(time, status) ~ education, fertility))
print(gray_test(Event(time, status) ~ education, fertility))

for (k in 1:2) {
  cat("\n---- cause", k, "----\n")
  fits <- list()
  for (b in c("gompertz", "mix_weibull", "mix_loglogistic")) for (l in c("PH", "PO", "GOR")) {
    f <- safe(directcif(Event(time, status) ~ education + age_c, fertility, baseline = b, link = l,
                        estimation = "weighted", cause = k))
    if (!is.null(f) && f@converged) fits[[paste(b, l)]] <- f
  }
  cat("weighted pseudo-log-likelihood by baseline and link:\n")
  print(data.frame(model = names(fits), loglik = vapply(fits, function(f) f@loglik, 0), npar = vapply(fits, function(f) length(coef(f)), 0)),
        row.names = FALSE)
  best <- fits[[which.max(vapply(fits, function(f) f@loglik - length(coef(f)), 0))]]
  cat("\nselected (largest pseudo-log-likelihood minus number of parameters):\n")
  print(best)
  cat("\nstandard errors with the censoring term (ipcw) and with the censoring distribution treated as known:\n")
  st <- se_tab(best); print(st, digits = 4)
  utils::write.csv(st, file.path(out_dir, paste0("fertility_cause", k, "_se.csv")), row.names = FALSE)
  nd <- data.frame(education = levels(fertility$education), age_c = 0)
  pr <- predict(best, nd, times = Inf, cause = k, se = TRUE)
  cat("\nlong-term incidence by education (age at the mean):\n")
  print(data.frame(education = nd$education, Finf = pr$cif[, 1], se = pr$se[, 1]), digits = 3)
  ff <- safe(fg(Event(time, status) ~ education + age_c, fertility, cause = k))
  if (!is.null(ff)) { cat("\nsemi-parametric Fine-Gray for comparison:\n"); print(ff) }
  pdf(file.path(out_dir, paste0("fertility_cause", k, ".pdf")), width = 7, height = 5)
  safe(cif_plot(best, cause = k, type = "unadjusted", reference = TRUE))
  dev.off()
}

# ------------------------------------------------------------------------------------------
# 2. flchain: three causes of death (circulatory, neoplasms, other) in a population cohort
#    aged 50 or more; everyone eventually dies, so the long-run risks add up to one.
# ------------------------------------------------------------------------------------------
cat("\n==== flchain ====\n")
data(flchain, package = "survival")
fl <- flchain
fl$time <- pmax(fl$futime, 1) / 365.25
ch <- as.character(fl$chapter)
fl$status <- ifelse(fl$death == 0, 0L, ifelse(ch %in% "Circulatory", 1L, ifelse(ch %in% "Neoplasms", 2L, 3L)))
fl$age10 <- (fl$age - 65) / 10
fl$male <- as.integer(fl$sex == "M")
fl$flc_high <- as.integer(fl$flc.grp >= 9)
cat("status: 0 alive, 1 circulatory, 2 neoplasms, 3 other\n"); print(table(fl$status))
print(aj(Event(time, status) ~ 1, fl))

form <- Event(time, status) ~ age10 + male + flc_high
nd <- expand.grid(age10 = c(0, 1), male = 0:1, flc_high = 0:1)
cons <- list()
for (b in c("gompertz", "mix_weibull", "logistic3")) for (g in c("weibull", "loglogistic", "lognormal")) {
  f <- safe(directcif(form, fl, baseline = b, estimation = "constrained", cause = 1, cond = g))
  if (!is.null(f) && f@converged) cons[[paste(b, g)]] <- f
}
cat("\nconstrained models, K = 3:\n")
print(data.frame(model = names(cons), loglik = vapply(cons, function(f) f@loglik, 0),
                 AIC = vapply(cons, function(f) AIC(f), 0)), row.names = FALSE)
if (length(cons)) {
  best <- cons[[which.min(vapply(cons, function(f) AIC(f), 0))]]
  print(summary(best))
  lr <- sapply(1:3, function(k) predict(best, nd, times = Inf, cause = k))
  colnames(lr) <- c("circulatory", "neoplasms", "other")
  cat("\nlong-run risks (constrained, best AIC) and their sum:\n")
  print(cbind(nd, round(lr, 3), sum = round(rowSums(lr), 6)))
  utils::write.csv(cbind(nd, lr), file.path(out_dir, "flchain_longrun_constrained.csv"), row.names = FALSE)
  cat("\nsensitivity: long-run risk of circulatory death at age 65, female, normal FLC, by model:\n")
  print(sapply(cons, function(f) predict(f, nd[1, ], times = Inf, cause = 1)))
}
fj <- safe(directcif(form, fl, baseline = "gompertz"))
if (!is.null(fj)) {
  cat("\njoint model without the constraint:\n"); print(fj)
  tot <- Reduce(`+`, lapply(1:3, function(k) predict(fj, fl, times = Inf, cause = k)))
  cat("share of subjects whose long-run risks add up to more than 1:", mean(tot > 1), "  largest sum:", max(tot), "\n")
}
fcs <- safe(csc(form, fl))
if (!is.null(fcs)) {
  t10 <- 10
  cat("\ncause-specific Cox, incidences at 10 years:\n")
  print(cbind(nd, sapply(1:3, function(k) predict(fcs, nd, times = t10, cause = k))))
  if (length(cons)) {
    cat("constrained model, incidences at 10 years:\n")
    print(cbind(nd, sapply(1:3, function(k) predict(best, nd, times = t10, cause = k))))
  }
}
if (length(cons)) {
  pdf(file.path(out_dir, "flchain_cif.pdf"), width = 9, height = 4)
  safe(cif_plot(best, cause = 1:3, type = "unadjusted", reference = TRUE, risk.table = FALSE))
  dev.off()
}
sink()
cat("written", file.path(out_dir, "realdata.txt"), "\n")
