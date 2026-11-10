# Open minor_update_fsusie.Rproj first, then source this script.
stopifnot(file.exists("config/study_config.R"))
source("R/study.R")
source("R/summarize.R")
cfg <- load_config()
cfg$run_id <- paste0(cfg$run_id, "_pilot")
cfg$n_rep <- 1L
run_study(cfg = cfg)
summarize_study(run_id = cfg$run_id)
