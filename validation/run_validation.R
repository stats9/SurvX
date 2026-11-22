#!/usr/bin/env Rscript
# SurvX validation suite: accuracy against established packages, speed, statistical properties
# (Monte Carlo) and stability. See README.md in this folder.
#
#   Rscript run_validation.R                          quick run of every part (about 20-40 minutes)
#   Rscript run_validation.R full cores=6 run=v2      the full study for the paper (several hours)
#   Rscript run_validation.R full cores=6 run=v2 parts=2     only some parts (1 accuracy, 2 speed,
#                                                     3 simulation, 4 stability), into the same run
#   Rscript run_validation.R run=v2 parts=report      only remake the report of run v2
#   Rscript run_validation.R full sim_R=500 sim_n=500 override any setting of make_config()
#
# The work is split into stages (accuracy: 2, speed: one per sample size and one for the number of
# covariates, simulation: one per sample size, stability: 1). Every finished stage is saved at once
# to validation_results/<run>/stages/, so a run can be stopped (Ctrl+C) and continued later with
# the same command: finished stages are skipped (redo=TRUE runs them again). Only the stage that
# was running when it stopped is lost. Without run=, a new run named by the date and time starts.
#
# After every session the report is made from all finished stages. Two files hold everything:
#   SurvX_validation_results.txt   all numeric results (summaries, every failed check, the stages
#                                  with their SurvX version and time, the log of every session)
#   SurvX_validation_figures.pdf   all figures
# The folder also keeps the full CSV files of every check (for the archive), tables/ (CSV and LaTeX
# for the paper) and validation_report.md (and .html with rmarkdown).

args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 error = function(e) getwd())
if (is.na(here) || !dir.exists(file.path(here, "R"))) here <- getwd()
WORKER_FILES <- normalizePath(file.path(here, "R", c("00_common.R", "03_simulation.R", "06_new_methods.R", "07_assess.R")))
for (f in c("00_common.R", "01_accuracy.R", "02_speed.R", "03_simulation.R", "04_stability.R", "05_report.R", "06_new_methods.R", "07_assess.R"))
  source(file.path(here, "R", f))

cfg <- make_config(args)
report_only <- identical(as.character(cfg$parts), "report")
if (report_only && is.null(cfg$run)) stop("parts=report needs the name of the run: run=<name>", call. = FALSE)
run_name <- if (is.null(cfg$run)) format(Sys.time(), "%Y%m%d-%H%M") else as.character(cfg$run)
out_dir <- file.path(getwd(), "validation_results", run_name)
stage_dir <- file.path(out_dir, "stages")
resumed <- dir.exists(stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
rec_init()

# ---- the settings of a run are fixed when it starts: a continued run must use the same ones -------
cfg_file <- file.path(out_dir, "config.csv")
keys <- c("mode", "seed", "sim_R", "sim_n", "speed_n", "speed_p", "speed_reps", "speed_max_seconds", "acc_reps")
now_cfg <- data.frame(setting = keys, value = vapply(keys, function(k) paste(cfg[[k]], collapse = ","), character(1)))
if (file.exists(cfg_file)) {
  old_cfg <- utils::read.csv(cfg_file, stringsAsFactors = FALSE, colClasses = "character")
  diff <- merge(old_cfg, now_cfg, by = "setting", suffixes = c("_run", "_now"))
  diff <- diff[diff$value_run != diff$value_now, , drop = FALSE]
  if (nrow(diff) && !report_only && !isTRUE(cfg$force))
    stop("run '", run_name, "' was started with other settings (", paste0(diff$setting, " = ", diff$value_run, collapse = "; "),
         "); use the same settings, another run=, or force=TRUE", call. = FALSE)
  for (k in keys) {                                   # the report describes the run as it was started
    v <- old_cfg$value[old_cfg$setting == k]
    if (length(v)) { num <- suppressWarnings(as.numeric(strsplit(v, ",")[[1]])); cfg[[k]] <- if (all(is.finite(num))) num else v }
  }
} else utils::write.csv(now_cfg, cfg_file, row.names = FALSE)

logmsg("SurvX ", format(packageVersion("SurvX")), " validation, run ", run_name, if (resumed) " (continued)" else " (new)",
       ", mode ", cfg$mode, ", parts ", paste(cfg$parts, collapse = ","), ", cores ", cfg$cores, "; output: ", out_dir)
logmsg("reference packages: ", paste(REF_PKGS, ifelse(vapply(REF_PKGS, has, logical(1)), "yes", "MISSING"), collapse = ", "))
cat("\n\n==== session ", format(Sys.time()), " ====\n", file = file.path(out_dir, "sessionInfo.txt"), append = TRUE)
cat(utils::capture.output(sessionInfo()), sep = "\n", file = file.path(out_dir, "sessionInfo.txt"), append = TRUE)

save_csv <- function(x, f, dir = out_dir) if (!is.null(x)) utils::write.csv(x, file.path(dir, f), row.names = FALSE)
read_csv <- function(f, dir = out_dir) { p <- file.path(dir, f); if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL }
bind <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(x) { for (c in setdiff(cols, names(x))) x[[c]] <- NA; x[cols] }))
}

# the log of every session goes to run.log.csv as it is written
log_done <- 0L
log_flush <- function() {
  lg <- rec_get("log")
  if (is.null(lg) || nrow(lg) <= log_done) return(invisible())
  new <- lg[(log_done + 1L):nrow(lg), , drop = FALSE]
  f <- file.path(out_dir, "run.log.csv")
  utils::write.table(new, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f), qmethod = "double")
  log_done <<- nrow(lg)
}

# one stage: skipped when its file exists (unless redo=TRUE); saved at once when finished
stage_rows <- list()
stage <- function(label, file, table, fun) {
  if (file.exists(file.path(stage_dir, file)) && !isTRUE(cfg$redo)) {
    logmsg("  stage ", .PROG$stage_k, " of ", .PROG$stage_n, ": ", label, ": finished earlier, skipped (", file, ")")
    return(invisible(read_csv(file, stage_dir)))
  }
  if (!is.null(table)) assign(table, list(), envir = .REC)
  t0 <- Sys.time()
  logmsg("  stage ", .PROG$stage_k, " of ", .PROG$stage_n, ": ", label)
  res <- tryCatch(fun(), error = function(e) { logmsg("  ", label, " stopped: ", conditionMessage(e)); NULL })
  ok <- !is.null(res) && NROW(res) > 0L
  if (ok) save_csv(res, file, stage_dir)
  row <- data.frame(stage = label, file = file, status = if (ok) "finished" else "failed", started = format(t0), finished = format(Sys.time()),
                    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), SurvX = format(packageVersion("SurvX")),
                    cores = cfg$cores, stringsAsFactors = FALSE)
  f <- file.path(out_dir, "stages.csv")
  utils::write.table(row, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f), qmethod = "double")
  log_flush()
  invisible(if (ok) res else NULL)
}

# the stages of every part, in order (label, file, recording table, function)
plan <- list(
  `1` = list(list("accuracy: established methods", "accuracy_main.csv", "accuracy", function() run_accuracy(cfg)),
             list("accuracy: methods added after 1.0.0", "accuracy_new.csv", "accuracy", function() run_accuracy_new(cfg))),
  `2` = c(lapply(cfg$speed_n, function(n) list(paste0("speed: n = ", n), sprintf("speed_n%07d.csv", as.integer(n)), "speed", function() {
            done <- bind(lapply(cfg$speed_n[cfg$speed_n < n], function(m) read_csv(sprintf("speed_n%07d.csv", as.integer(m)), stage_dir)))
            run_speed_stage(cfg, "n", n, speed_too_slow(done, cfg$speed_max_seconds))
          })),
          list(list("speed: number of covariates", "speed_p.csv", "speed", function() run_speed_stage(cfg, "p")))),
  `3` = lapply(cfg$sim_n, function(n) list(paste0("simulation: n = ", n), sprintf("simulation_n%06d.csv", as.integer(n)), NULL,
                                          function() run_simulation_stage(cfg, n))),
  `4` = list(list("stability and invariances", "stability.csv", "stability", function() run_stability(cfg))))

t_all <- proc.time()[["elapsed"]]
if (!report_only) {
  sel <- intersect(c("1", "2", "3", "4"), as.character(cfg$parts))
  n_stages <- sum(lengths(plan[sel])); k_stage <- 0L
  for (part in sel) {
    logmsg(c(`1` = "Part 1: accuracy against reference implementations", `2` = "Part 2: speed",
             `3` = paste0("Part 3: Monte Carlo study (", cfg$sim_R, " replications per scenario, n = ", paste(cfg$sim_n, collapse = ", "), ")"),
             `4` = "Part 4: stability and invariances")[[part]])
    for (st in plan[[part]]) {
      k_stage <- k_stage + 1L
      progress_stage(k_stage, n_stages, st[[1]])
      stage(st[[1]], st[[2]], st[[3]], st[[4]])
    }
  }
}

# ---- the tables of the report, from every finished stage -----------------------------------------
files_of <- function(part) vapply(plan[[part]], `[[`, character(1), 2L)
missing <- unlist(lapply(names(plan), function(p) { f <- files_of(p); lab <- vapply(plan[[p]], `[[`, character(1), 1L)
  lab[!file.exists(file.path(stage_dir, f))] }))
save_csv(bind(lapply(files_of("1"), read_csv, dir = stage_dir)), "accuracy.csv")
save_csv(bind(lapply(files_of("2"), read_csv, dir = stage_dir)), "speed.csv")
raw <- bind(lapply(files_of("3"), read_csv, dir = stage_dir))
if (!is.null(raw)) {
  save_csv(raw, "simulation_raw.csv")
  save_csv(tryCatch(summarise_simulation(raw), error = function(e) { logmsg("simulation summary failed: ", conditionMessage(e)); NULL }),
           "simulation_summary.csv")
}
save_csv(read_csv("stability.csv", stage_dir), "stability.csv")
logmsg(if (length(missing)) paste0("stages not finished yet: ", paste(missing, collapse = "; ")) else "all stages finished")
log_flush()
verdict <- tryCatch(make_report(cfg, out_dir, missing), error = function(e) { logmsg("report failed: ", conditionMessage(e)); "report failed" })
logmsg("session finished in ", round((proc.time()[["elapsed"]] - t_all) / 60, 1), " minutes")
log_flush()
cat("\nVERDICT\n", if (length(verdict)) paste("-", verdict, collapse = "\n") else "- all checks passed", "\n",
    if (length(missing)) paste0("\nNot finished yet (continue with the same command, run=", run_name, "):\n",
                                paste("  -", missing, collapse = "\n"), "\n") else "\nAll stages finished.\n",
    "\nFor review, send these two files:\n  ", file.path(out_dir, "SurvX_validation_results.txt"),
    "\n  ", file.path(out_dir, "SurvX_validation_figures.pdf"), "\n", sep = "")
