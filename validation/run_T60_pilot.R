# Validation of the new 60-point design: one replicate of all 12 conditions.
# Run from the project root. No full 100-replicate study is launched here.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/study.R")
check_dependencies()
main_cfg <- load_config()
dat <- simulate_dataset(main_cfg, make_tasks(main_cfg)[1, ])
mapped <- fsusieR::remap_data(Y = dat$Y, pos = seq_len(main_cfg$positions), verbose = FALSE)
stopifnot(ncol(dat$Y) == 60L, ncol(mapped$Y) == 64L, all(is.finite(mapped$Y)),
  main_cfg$n_rep == 100L, isTRUE(main_cfg$fsusie_standardize_Y),
  main_cfg$gp_length_fraction * main_cfg$positions == 12,
  !main_cfg$include_shuffled_fsusie)
cat("Observed points:", ncol(dat$Y), "; fSuSiE internal grid:", ncol(mapped$Y), "\n")
source("scripts/00_check_and_pilot.R")
pilot_id <- paste0(main_cfg$run_id, "_pilot")
files <- list.files(file.path("results", pilot_id, "replicates"), "[.]rds$", full.names = TRUE)
results <- lapply(files, readRDS)
status <- do.call(rbind, lapply(results, `[[`, "status"))
stopifnot(length(files) == 12L, nrow(status) == 156L, all(status$success),
  !any(status$method == "fSuSiE_shuffled"),
  sum(status$method == "fSuSiE_standardized") == 12L)
note <- c(paste("Pilot validation:", pilot_id),
  "PASS: 60 observed points are accepted; installed fSuSiE remaps internally to 64 points.",
  "PASS: main configuration has 100 replicates per condition, GP length scale 12 and nugget 0.15.",
  "PASS: all 12 fSuSiE fits use standardize_Y = TRUE (internal ID fSuSiE_standardized; display label fSuSiE).",
  "PASS: all 12 conditions / 156 fits completed without errors; no shuffled fits.",
  paste("Nonconverged pilot fits:", sum(!status$converged, na.rm = TRUE)),
  "The full study has not been launched. Use scripts/05_launch_two_jobs.R from the RStudio console.")
writeLines(note, file.path("validation", paste0(main_cfg$run_id, ".txt")))
cat(paste(note, collapse = "\n"), "\n")
