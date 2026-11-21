# Paper 2: tables of the simulation study.
#   Rscript paper2/sim/summarize.R quick      (or full)
# Writes paper2/results/tables_<mode>.txt and the CSV files behind them.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "quick"
res <- readRDS(file.path("paper2", "results", paste0("sim_", mode, ".rds")))
z <- qnorm(0.975)

perf <- function(d, truth, prob) {
  ok <- is.finite(d$est)
  d <- d[ok, , drop = FALSE]
  th <- truth[[d$estimand[1]]]
  hasse <- is.finite(d$se)
  cover <- if (any(hasse)) {
    if (prob) {
      g <- log(-log1p(-d$est)); sg <- d$se / ((1 - d$est) * (-log1p(-d$est)))
      lo <- -expm1(-exp(g - z * sg)); hi <- -expm1(-exp(g + z * sg))
    } else { lo <- d$est - z * d$se; hi <- d$est + z * d$se }
    mean(lo[hasse] <= th & th <= hi[hasse], na.rm = TRUE)
  } else NA_real_
  sd <- stats::sd(d$est)
  data.frame(reps = nrow(d), true = th, bias = mean(d$est) - th, emp_sd = sd,
             se_sd = if (any(hasse)) mean(d$se[hasse]) / sd else NA_real_, rmse = sqrt(mean((d$est - th)^2)),
             coverage = cover)
}

tab <- function(D, scen, truth) {
  out <- list()
  for (s in seq_len(nrow(scen))) {
    Ds <- D[D$scenario == s, , drop = FALSE]
    for (e in unique(Ds$estimator)) for (m in unique(Ds$estimand[Ds$estimator == e])) {
      d <- Ds[Ds$estimator == e & Ds$estimand == m, , drop = FALSE]
      if (!nrow(d) || all(!is.finite(d$est))) next
      out[[length(out) + 1L]] <- cbind(scen[s, , drop = FALSE], estimator = e, estimand = m,
                                       perf(d, truth[[s]], !startsWith(m, "b_")), row.names = NULL)
    }
  }
  do.call(rbind, out)
}

TA <- tab(res$A, res$scen_A, res$truth_A)
TB <- tab(res$B, res$scen_B, res$truth_B)
coh <- stats::aggregate(cbind(excess_inf, excess_end, max_inf) ~ scenario + estimator, data = res$B_coherence, FUN = mean)
coh <- merge(cbind(scenario = seq_len(nrow(res$scen_B)), res$scen_B), coh)

dir.create(file.path("paper2", "results"), showWarnings = FALSE)
utils::write.csv(TA, file.path("paper2", "results", paste0("tableA_", mode, ".csv")), row.names = FALSE)
utils::write.csv(TB, file.path("paper2", "results", paste0("tableB_", mode, ".csv")), row.names = FALSE)
utils::write.csv(coh, file.path("paper2", "results", paste0("tableB_coherence_", mode, ".csv")), row.names = FALSE)
utils::write.csv(res$C, file.path("paper2", "results", paste0("tableC_", mode, ".csv")), row.names = FALSE)

sink(file.path("paper2", "results", paste0("tables_", mode, ".txt")))
cat("SurvX", res$version, "| mode", res$mode, "| replications per scenario", res$R, "\n\n")
op <- options(width = 200, digits = 3)
cat("Study A: weighted estimator, variance with and without the censoring term\n")
print(TA, row.names = FALSE)
cat("\nStudy B: three causes\n")
print(TB, row.names = FALSE)
cat("\nStudy B: share of the sample with predicted incidences adding up to more than 1\n")
print(coh, row.names = FALSE)
cat("\nStudy C: computing time (seconds)\n")
print(res$C, row.names = FALSE)
options(op)
sink()
cat("written", file.path("paper2", "results", paste0("tables_", mode, ".txt")), "\n")
