# OPTIONAL: launch both scripts as separate RStudio Background Jobs.
# This file does not run automatically; source it from the project console.
stopifnot(file.exists("config/study_config.R"))
if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable())
  stop("Run this launcher inside RStudio, or use scripts 01 and 02 in separate R sessions.")
root <- normalizePath(getwd(), winslash = "/")
for (name in c("01_run_worker_1.R", "02_run_worker_2.R")) {
  rstudioapi::jobRunScript(path = file.path(root, "scripts", name),
    name = sub("[.]R$", "", name), workingDir = root,
    importEnv = FALSE, exportEnv = "")
}
