# Run after both workers finish. Also safe for an explicitly labeled preview.
stopifnot(file.exists("config/study_config.R"))
source("R/study.R")
source("R/summarize.R")
summarize_study()
