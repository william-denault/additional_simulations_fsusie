# Run in the second RStudio session (or as a second Background Job).
stopifnot(file.exists("config/study_config.R"))
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/study.R")
run_study(worker_id = 2L, n_workers = 2L)
