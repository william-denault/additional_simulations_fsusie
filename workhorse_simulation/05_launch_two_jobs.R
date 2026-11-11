# SOURCE THIS FILE IN THE MAIN RSTUDIO CONSOLE, not as a Background Job.
# It creates two Background Jobs; scripts 01 and 02 can also be launched manually.
# Source this whole file, or use RStudio's "Source as Background Job".
.project_root = local({
  files = vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  files = files[nzchar(files)]
  args = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
  # RStudio Background Jobs may expose neither ofile nor --file. The launcher
  # supplies workingDir, so the current directory is a valid fallback.
  candidates = unique(c(dirname(rev(files)),dirname(args),getwd()))
  is_project = vapply(candidates,function(path){
    all(file.exists(file.path(path,c("workhorse_simulation.Rproj","config.R","R/setup.R"))))
  },logical(1))
  if (!any(is_project))
    stop("Cannot locate workhorse_simulation. Set the Background Job working directory to the project folder, or launch it with 05_launch_two_jobs.R.")
  normalizePath(candidates[which(is_project)[1]],winslash="/",mustWork=TRUE)
})
source(file.path(.project_root,"R","setup.R"))
setup(.project_root)
if (!requireNamespace("rstudioapi",quietly=TRUE) || !rstudioapi::isAvailable())
  stop("Source 05_launch_two_jobs.R in the main RStudio console, or run scripts 01 and 02 as separate Background Jobs.")
cfg = read_config(.project_root)
for (worker in 1:2){
  script = sprintf("%02d_run_worker_%d.R",worker,worker)
  rstudioapi::jobRunScript(path=file.path(.project_root,script),
    name=paste(cfg$run_id,paste0("[",cfg$n_rep," reps]"),"worker",worker),workingDir=.project_root,
    importEnv=FALSE,exportEnv="")
}
