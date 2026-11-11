# Optional: start/resume the WORKHORSE study as two RStudio Background Jobs.
# Source this whole script from the MAIN RStudio console, not as a Background Job.
local({
  files = vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  files = files[nzchar(files)]
  args = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
  candidates = unique(c(dirname(rev(files)), dirname(args), getwd(), dirname(getwd())))
  valid = vapply(candidates, function(path) all(file.exists(file.path(path,
    c("config/study_config.R", "workhorse_simulation/05_launch_two_jobs.R")))), logical(1))
  if (!any(valid)) stop("Source the whole script or open minor_update_fsusie.Rproj first.")
  root = normalizePath(candidates[which(valid)[1]], winslash="/", mustWork=TRUE)
  source(file.path(root, "workhorse_simulation", "05_launch_two_jobs.R"), local=TRUE)
})
