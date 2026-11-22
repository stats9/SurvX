#!/usr/bin/env Rscript
# SurvX validation suite: accuracy against established packages, speed, statistical properties
# (Monte Carlo) and stability. See README.md in this folder.
#
#   Rscript run_validation.R                 quick run of every part (about 20-40 minutes)
#   Rscript run_validation.R full cores=6    the full study for the paper (several hours)
#   Rscript run_validation.R parts=1,4       only some parts (1 accuracy, 2 speed, 3 simulation, 4 stability)
#   Rscript run_validation.R full sim_R=500 sim_n=500   override any setting of make_config()
#
# Results go to validation_results/<date-time>/. Two files hold everything needed for review:
#   SurvX_validation_results.txt   all numeric results (summaries, every failed check, log, session)
#   SurvX_validation_figures.pdf   all figures
# The folder also keeps the full CSV files of every check (for the archive), tables/ (CSV and LaTeX
# for the paper) and validation_report.md (and .html with rmarkdown).

args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 error = function(e) getwd())
if (is.na(here) || !dir.exists(file.path(here, "R"))) here <- getwd()
WORKER_FILES <- normalizePath(file.path(here, "R", c("00_common.R", "03_simulation.R", "06_new_methods.R")))
for (f in c("00_common.R", "01_accuracy.R", "02_speed.R", "03_simulation.R", "04_stability.R", "05_report.R", "06_new_methods.R"))
  source(file.path(here, "R", f))

cfg <- make_config(args)
out_dir <- file.path(getwd(), "validation_results", format(Sys.time(), "%Y%m%d-%H%M"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rec_init()
logmsg("SurvX ", format(packageVersion("SurvX")), " validation, mode ", cfg$mode, ", parts ", paste(cfg$parts, collapse = ","),
       ", cores ", cfg$cores, "; output: ", out_dir)
logmsg("reference packages: ", paste(REF_PKGS, ifelse(vapply(REF_PKGS, has, logical(1)), "yes", "MISSING"), collapse = ", "))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
save_csv <- function(x, f) if (!is.null(x)) utils::write.csv(x, file.path(out_dir, f), row.names = FALSE)

t_all <- proc.time()[["elapsed"]]
if (1 %in% cfg$parts) save_csv(tryCatch(run_accuracy(cfg), error = function(e) { logmsg("accuracy stopped: ", conditionMessage(e)); rec_get("accuracy") }), "accuracy.csv")
if (2 %in% cfg$parts) save_csv(tryCatch(run_speed(cfg), error = function(e) { logmsg("speed stopped: ", conditionMessage(e)); rec_get("speed") }), "speed.csv")
if (3 %in% cfg$parts) {
  raw <- tryCatch(run_simulation(cfg), error = function(e) { logmsg("simulation stopped: ", conditionMessage(e)); NULL })
  if (!is.null(raw)) {
    save_csv(raw, "simulation_raw.csv")
    save_csv(tryCatch(summarise_simulation(raw), error = function(e) { logmsg("simulation summary failed: ", conditionMessage(e)); NULL }),
             "simulation_summary.csv")
  }
}
if (4 %in% cfg$parts) save_csv(tryCatch(run_stability(cfg), error = function(e) { logmsg("stability stopped: ", conditionMessage(e)); rec_get("stability") }), "stability.csv")
verdict <- tryCatch(make_report(cfg, out_dir), error = function(e) { logmsg("report failed: ", conditionMessage(e)); "report failed" })
save_csv(rec_get("log"), "run.log.csv")
logmsg("finished in ", round((proc.time()[["elapsed"]] - t_all) / 60, 1), " minutes")
cat("\nVERDICT\n", if (length(verdict)) paste("-", verdict, collapse = "\n") else "- all checks passed", "\n",
    "\nFor review, send these two files:\n  ", file.path(out_dir, "SurvX_validation_results.txt"),
    "\n  ", file.path(out_dir, "SurvX_validation_figures.pdf"), "\n", sep = "")
