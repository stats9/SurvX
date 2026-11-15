#!/usr/bin/env Rscript
# =====================================================================================================
#  Simulation study 2: the published models and the extensions implemented in SurvX
#
#  Everything of sim_published_models.R, plus
#    Weighted PH / PO / GOR    cause 1 alone; competing events through Kaplan-Meier censoring weights
#                              as in Fine & Gray (1999), Gompertz baseline, sandwich variance
#                              (directcif(estimation = "weighted"))
#    Weighted (AIC)            the weighted fit with the smallest pseudo-AIC
#    Constrained PH / PO, G by AIC
#                              Shi, Cheng & Jeong (2013) with G chosen by AIC among the modified
#                              logistic, Weibull, log-logistic and log-normal forms
#    Direct spline PH / PO     joint model with a restricted cubic spline baseline (3 df)
#    SurvX (AIC)               the full-likelihood fit (joint Gompertz or spline, constrained) with the
#                              smallest AIC: an automatic choice an analyst can make from the data
#
#  Run (all settings optional; same seeds as study 1, so the common models give identical results):
#    Rscript sim_survx_models.R R=1000 cores=8 sizes=100,250,500,1000 scenarios=S1,S2,S3,S4,S5,S6,S7,S8
#  Trial:  Rscript sim_survx_models.R R=20 sizes=250
#  This study fits 17 models per data set and takes about three times as long as study 1.
# =====================================================================================================

# ---- SETTINGS ---------------------------------------------------------------------------------------
SETTINGS <- list(
  R         = 1000L,
  sizes     = c(100L, 250L, 500L, 1000L),
  scenarios = paste0("S", 1:8),
  cores     = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  n_test    = 10000L,
  seed      = 20261115L,
  out_dir   = "results"
)
# ------------------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L || !kv[1] %in% names(SETTINGS)) stop("unknown argument: ", a)
  SETTINGS[[kv[1]]] <- switch(kv[1], scenarios = strsplit(kv[2], ",")[[1]], out_dir = kv[2],
                              as.integer(strsplit(kv[2], ",")[[1]]))
}

core_file <- file.path(local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}), "sim_core.R")
if (!file.exists(core_file)) stop("sim_core.R not found next to this script; run the script from its own folder")
source(core_file)

G_all <- c("logistic", "weibull", "loglogistic", "lognormal")
specs_survx <- list(
  # published (identical to study 1)
  list(label = "Direct PH",          type = "joint", link = "PH",  group = "direct (published)"),
  list(label = "Direct PO",          type = "joint", link = "PO",  group = "direct (published)"),
  list(label = "Direct GOR",         type = "joint", link = "GOR", group = "direct (published)"),
  list(label = "Direct (AIC)",       type = "select", from = c("Direct PH", "Direct PO", "Direct GOR"), group = "direct (published)"),
  list(label = "Constrained PH",     type = "constrained", link = "PH", cond = "logistic", group = "direct (published)"),
  list(label = "Constrained PO",     type = "constrained", link = "PO", cond = "logistic", group = "direct (published)"),
  # SurvX extensions
  list(label = "Weighted PH",        type = "weighted", link = "PH",  group = "direct (SurvX)"),
  list(label = "Weighted PO",        type = "weighted", link = "PO",  group = "direct (SurvX)"),
  list(label = "Weighted GOR",       type = "weighted", link = "GOR", group = "direct (SurvX)"),
  list(label = "Weighted (AIC)",     type = "select", from = c("Weighted PH", "Weighted PO", "Weighted GOR"), group = "direct (SurvX)"),
  list(label = "Constrained PH, G by AIC", type = "constrained", link = "PH", cond = G_all, group = "direct (SurvX)"),
  list(label = "Constrained PO, G by AIC", type = "constrained", link = "PO", cond = G_all, group = "direct (SurvX)"),
  list(label = "Direct spline PH",   type = "joint", link = "PH", baseline = "spline", df = 3L, group = "direct (SurvX)"),
  list(label = "Direct spline PO",   type = "joint", link = "PO", baseline = "spline", df = 3L, group = "direct (SurvX)"),
  list(label = "SurvX (AIC)",        type = "select",
       from = c("Direct PH", "Direct PO", "Direct GOR", "Direct spline PH", "Direct spline PO",
                "Constrained PH, G by AIC", "Constrained PO, G by AIC"), group = "direct (SurvX)"),
  # references
  list(label = "Fine-Gray",          type = "fg",  group = "reference"),
  list(label = "Cause-specific Cox", type = "csc", group = "reference")
)

cat("SurvX", format(utils::packageVersion("SurvX")), "| study: SurvX models |", SETTINGS$R, "replications |",
    "sizes", paste(SETTINGS$sizes, collapse = ","), "| scenarios", paste(SETTINGS$scenarios, collapse = ","),
    "|", SETTINGS$cores, "cores\n\n")
raw <- run_study("survx", specs_survx, scenario_ids = SETTINGS$scenarios, sizes = SETTINGS$sizes, R = SETTINGS$R,
                 n_test = SETTINGS$n_test, seed = SETTINGS$seed, cores = SETTINGS$cores, out_dir = SETTINGS$out_dir,
                 core_file = core_file)
summ <- summarise_study(raw)
utils::write.csv(summ, file.path(SETTINGS$out_dir, "summary_survx.csv"), row.names = FALSE)
utils::write.csv(summarise_design(raw), file.path(SETTINGS$out_dir, "design_survx.csv"), row.names = FALSE)
saveRDS(list(specs = specs_survx, settings = SETTINGS, scenarios = scenarios()[SETTINGS$scenarios],
             session = utils::sessionInfo()), file.path(SETTINGS$out_dir, "setup_survx.rds"))
cat("written:", file.path(SETTINGS$out_dir, c("summary_survx.csv", "design_survx.csv")), "\n")
cat("next: Rscript sim_figures.R\n")
