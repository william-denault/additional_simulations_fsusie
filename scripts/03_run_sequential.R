# Alternative to the TWO worker scripts. Do not start both workflows together.
stopifnot(file.exists("config/study_config.R"))
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/study.R")
run_study()
