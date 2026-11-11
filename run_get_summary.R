# Regenerate BOTH completed studies from saved checkpoints; no model fitting.
# Source the whole script, or run as a Background Job from the repository root.
local({
  files = vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  files = files[nzchar(files)]
  args = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
  candidates = unique(c(dirname(rev(files)), dirname(args), getwd(), dirname(getwd())))
  valid = vapply(candidates, function(path) all(file.exists(file.path(path,
    c("config/study_config.R", "scripts/04_summarize.R",
      "workhorse_simulation/04_summarize_and_plot.R")))), logical(1))
  if (!any(valid)) stop("Source the whole script, or use the repository root as the Background Job working directory.")
  root = normalizePath(candidates[which(valid)[1]], winslash="/", mustWork=TRUE)
  previous = setwd(root)
  on.exit(setwd(previous), add=TRUE)

  message("Regenerating covariance benchmark: 100 replicates per condition.")
  source(file.path(root, "scripts", "04_summarize.R"), local=TRUE)
  message("Regenerating workhorse: 250 replicates per condition.")
  source(file.path(root, "workhorse_simulation", "04_summarize_and_plot.R"), local=TRUE)
  message("Done. Both studies' tables and plots were regenerated from saved results.")
})
