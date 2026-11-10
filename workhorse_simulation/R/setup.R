setup = function(root){
  Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1)
  if (.Platform$OS.type=="windows"){
    v = paste(R.version$major, strsplit(R.version$minor,".",fixed=TRUE)[[1]][1],sep=".")
    libs = c(file.path(Sys.getenv("USERPROFILE"),"AppData","Local","R","win-library",v),
             file.path(Sys.getenv("LOCALAPPDATA"),"R","win-library",v))
    .libPaths(unique(c(.libPaths(), libs[dir.exists(libs)])))
  }
  packages = c("susieR","fsusieR","mvsusieR","mashr")
  absent = packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(absent)) stop("Install these packages first: ",paste(absent,collapse=", "))
  for (file in c("get_sim.R","run_methods.R","study.R","plots.R"))
    sys.source(file.path(root,"R",file), envir=.GlobalEnv)
  invisible(root)
}
