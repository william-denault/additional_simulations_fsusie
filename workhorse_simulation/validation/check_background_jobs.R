# Run from the project directory. No actual simulations/jobs are started.
root = normalizePath(getwd(),winslash="/")
files = list.files(root,pattern="^0[0-5]_.*\\.R$",full.names=TRUE)
stopifnot(length(files)==6)

job_environment = function(){
  e = new.env(parent=baseenv())
  # RStudio can evaluate expressions without source() metadata or --file.
  e$sys.frames = function() list()
  e$commandArgs = function(...) c("Rterm","--no-save")
  e
}

for (file in files){
  e = job_environment()
  eval(parse(file)[[1]],envir=e)
  stopifnot(identical(e$.project_root,root))
  # Ignore a wrapper's --file path if it is not the workhorse project.
  e = job_environment()
  e$commandArgs = function(...) c("Rterm",paste0("--file=",file.path(root,"validation","job_wrapper.R")))
  eval(parse(file)[[1]],envir=e)
  stopifnot(identical(e$.project_root,root))
}

for (worker in 1:2){
  e = job_environment()
  e$source = function(...) invisible(NULL)
  e$setup = function(...) invisible(NULL)
  e$run_study = function(root,worker,n_workers) list(root=root,worker=worker,n_workers=n_workers)
  result = eval(parse(file.path(root,sprintf("%02d_run_worker_%d.R",worker,worker))),envir=e)
  stopifnot(result$root==root,result$worker==worker,result$n_workers==2L)
}

# Sourcing by absolute path must still work from an unrelated working folder.
for (file in files){
  e = job_environment()
  e$sys.frames = local({
    path = file
    function() list(list(ofile=path))
  })
  e$getwd = function() file.path(root,"validation")
  eval(parse(file)[[1]],envir=e)
  stopifnot(identical(e$.project_root,root))
}
cat("PASS: all six entry scripts locate the project without source-file metadata,\n",
    "ignore unrelated wrapper paths, preserve absolute-path sourcing, and both workers dispatch correctly.\n")
