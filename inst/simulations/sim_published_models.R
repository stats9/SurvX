#!/usr/bin/env Rscript
# =====================================================================================================
#  Simulation study 1: published models only
#
#  Direct parametric models of the cumulative incidence function as published, against Fine-Gray
#  and cause-specific Cox models:
#    Direct PH / PO / GOR   joint likelihood of all causes, Gompertz baseline
#                           (Jeong & Fine 2006, 2007)
#    Direct (AIC)           the one of the three with the smallest AIC in each data set
#    Constrained PH / PO    two causes whose long-run risks add up to one, cause 2 through the
#                           modified logistic G (Shi, Cheng & Jeong 2013)
#    Fine-Gray              Fine & Gray (1999)
#    Cause-specific Cox     one Cox model per cause, combined into the CIF (Benichou & Gail 1990)
#
#  Run from a terminal in this directory (all settings optional):
#    Rscript sim_published_models.R R=1000 cores=8 sizes=100,250,500,1000 scenarios=S1,S2,S3,S4,S5,S6,S7,S8
#  or in RStudio: set the working directory to this folder, edit SETTINGS and source the file.
#
#  A short trial first (a few minutes):   Rscript sim_published_models.R R=20 sizes=250
#  (cells are stored with the number of replications in their name, so the trial does not mix
#  with the real run).
#
#  Output in results/: cells/*.rds (one file per scenario x sample size x R; finished cells are skipped
#  when the script is restarted), raw_published.csv.gz, summary_published.csv, design_published.csv.
#  Figures and tables: Rscript sim_figures.R
# =====================================================================================================

# ---- SETTINGS ---------------------------------------------------------------------------------------
SETTINGS <- list(
  R         = 1000L,                                   # replications per cell
  sizes     = c(100L, 250L, 500L, 1000L),              # training sample sizes
  scenarios = paste0("S", 1:8),                        # see sim_core.R for the full description
  cores     = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),  # parallel workers
  n_test    = 10000L,                                  # test subjects per scenario
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

specs_published <- list(
  list(label = "Direct PH",          type = "joint", link = "PH",  group = "direct (published)"),
  list(label = "Direct PO",          type = "joint", link = "PO",  group = "direct (published)"),
  list(label = "Direct GOR",         type = "joint", link = "GOR", group = "direct (published)"),
  list(label = "Direct (AIC)",       type = "select", from = c("Direct PH", "Direct PO", "Direct GOR"), group = "direct (published)"),
  list(label = "Constrained PH",     type = "constrained", link = "PH", cond = "logistic", group = "direct (published)"),
  list(label = "Constrained PO",     type = "constrained", link = "PO", cond = "logistic", group = "direct (published)"),
  list(label = "Fine-Gray",          type = "fg",  group = "reference"),
  list(label = "Cause-specific Cox", type = "csc", group = "reference")
)

cat("SurvX", format(utils::packageVersion("SurvX")), "| study: published models |", SETTINGS$R, "replications |",
    "sizes", paste(SETTINGS$sizes, collapse = ","), "| scenarios", paste(SETTINGS$scenarios, collapse = ","),
    "|", SETTINGS$cores, "cores\n\n")
raw <- run_study("published", specs_published, scenario_ids = SETTINGS$scenarios, sizes = SETTINGS$sizes, R = SETTINGS$R,
                 n_test = SETTINGS$n_test, seed = SETTINGS$seed, cores = SETTINGS$cores, out_dir = SETTINGS$out_dir,
                 core_file = core_file)
summ <- summarise_study(raw)
utils::write.csv(summ, file.path(SETTINGS$out_dir, "summary_published.csv"), row.names = FALSE)
utils::write.csv(summarise_design(raw), file.path(SETTINGS$out_dir, "design_published.csv"), row.names = FALSE)
saveRDS(list(specs = specs_published, settings = SETTINGS, scenarios = scenarios()[SETTINGS$scenarios],
             session = utils::sessionInfo()), file.path(SETTINGS$out_dir, "setup_published.rds"))
cat("written:", file.path(SETTINGS$out_dir, c("summary_published.csv", "design_published.csv")), "\n")
cat("next: Rscript sim_figures.R\n")
