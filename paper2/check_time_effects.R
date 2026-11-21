# Quick check of time_effects() before CI. Run from the repository root after installing dev:
#   R CMD INSTALL .
#   Rscript paper2/check_time_effects.R
# Writes check_time_effects.txt/.pdf/_report.html and check_adjusted_risk.txt/.pdf/_report.html.
library(SurvX)
sink("check_time_effects.txt", split = TRUE)
pdf("check_time_effects.pdf", width = 8, height = 6)
cat("SurvX", as.character(packageVersion("SurvX")), "\n")
step <- function(label, expr) {
  cat("\n=====", label, "=====\n")
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(expr, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  cat("(", round(proc.time()[["elapsed"]] - t0, 2), "s )\n")
  invisible(r)
}

# 1. Cox, veteran: karno is the classic time-varying effect
vet <- survival::veteran
fv <- cox(Event(time, status) ~ trt + celltype + karno, vet)
te <- step("cox veteran", { x <- time_effects(fv, times = c(30, 90, 180, 365)); print(x); x })
step("summary", summary(te))
step("plot", plot(te))
step("as.data.frame forms", print(as.data.frame(te, what = "forms")))
step("predict", print(predict(te, data.frame(trt = 1, celltype = "squamous", karno = c(40, 80)), times = c(30, 90, 180))))
step("fix piecewise at 90 days", print(time_effects(fv, fix = list(karno = "piecewise"), breaks = 90)))
step("own form sqrt", print(time_effects(fv, forms = list("constant", "log", sqrt = function(t) sqrt(t)))))
step("bootstrap 50", { b <- time_effects(fv, bootstrap = 50); print(b); report(b, file = "check_time_effects_report.html"); b })

# 2. agreement with survival::coxph(tt = ) for the log form
step("log form against coxph(tt)", {
  x <- time_effects(fv, forms = c("constant", "log"), fix = list(karno = "log"))
  et <- sort(unique(vet$time[vet$status == 1])); tm <- median(et)
  r <- survival::coxph(survival::Surv(time, status) ~ trt + celltype + karno + tt(karno), data = vet,
                       tt = function(x, t, ...) x * log(t / tm))
  print(cbind(SurvX = x@units[[1]]$final$coefficients, survival = coef(r)))
  cat("max abs difference:", max(abs(x@units[[1]]$final$coefficients - coef(r))), "\n")
})

# 3. Fine-Gray and cause-specific Cox on fertility
data(fertility)
step("fg fertility", print(time_effects(fg(Event(time, status) ~ education, fertility, cause = 1))))
step("csc fertility", print(time_effects(csc(Event(time, status) ~ education, fertility))))

# 4. flchain (larger data): Cox for death, age and sex
fl <- survival::flchain
fl$time <- pmax(fl$futime, 1) / 365.25
fl$male <- as.integer(fl$sex == "M"); fl$age10 <- (fl$age - 65) / 10
step("cox flchain", print(time_effects(cox(Event(time, death) ~ age10 + male, fl))))

dev.off()
sink()

# 5. adjusted_risk()
sink("check_adjusted_risk.txt", split = TRUE)
pdf("check_adjusted_risk.pdf", width = 9, height = 4.5)
lung <- survival::lung
lung$status <- lung$status - 1; lung$sex <- factor(lung$sex, 1:2, c("male", "female")); lung$years <- lung$time / 365.25
fl <- cox(Event(years, status) ~ sex + age + ph.ecog, data = lung)
ar <- step("lung: both methods, B = 0", { a <- adjusted_risk(fl, lung, "sex", times = c(0.5, 1), method = "both", B = 0); print(a); plot(a); a })
step("lung: bootstrap 200", { a <- adjusted_risk(fl, lung, "sex", times = c(0.5, 1), method = "both", B = 200); summary(a)
                              report(a, file = "check_adjusted_risk_report.html"); a })
step("fertility fg", print(adjusted_risk(fg(Event(time, status) ~ education, data = fertility, cause = 1), fertility, "education",
                                         times = c(1, 2, 5), method = "both", B = 100)))
step("fertility csc", print(adjusted_risk(csc(Event(time, status) ~ education, data = fertility), fertility, "education",
                                          times = c(1, 2, 5), B = 0)))
step("time_effects fit", print(adjusted_risk(time_effects(fv), vet, "trt", times = c(30, 90, 180), B = 0)))
dev.off()
sink()
