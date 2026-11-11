# Check extending an existing run without changing or refitting old data.
# Temporary synthetic checkpoint fixtures are removed after this check.
root = normalizePath(getwd(),winslash="/")
source(file.path(root,"R","setup.R")); setup(root)
check_extension = function(){
  cfg = read_config(root)
  cfg$n_rep = 250L
  old_cfg = cfg; old_cfg$n_rep = 50L
  grid = study_grid(cfg); old_grid = study_grid(old_cfg)
  sig = study_signature(root,cfg)
  old_sig = study_signature(root,old_cfg)
  old_sig$config$n_rep = 50L # Reproduce the legacy signature format.
  stopifnot(compatible_signatures(old_sig,sig),nrow(grid)==2000,nrow(old_grid)==400,
    all(vapply(seq_len(nrow(old_grid)),function(i) same_task(old_grid[i,,drop=FALSE],grid[i,,drop=FALSE]),TRUE)),
    length(intersect(unique(old_grid$seed),unique(grid$seed[grid$replicate>50])))==0)
  for (field in c("seed","n_snps","L_fit","gp_length_scale")){
    bad = sig; bad$config[[field]] = bad$config[[field]]+1
    stopifnot(!compatible_signatures(old_sig,bad))
  }
  bad = sig; bad$code[1] = "changed fitting code"
  stopifnot(!compatible_signatures(old_sig,bad))
  bad = sig; bad$packages[1] = "changed package version"
  stopifnot(!compatible_signatures(old_sig,bad))
  for (worker in 1:2){
    tasks = grid[(grid$replicate-1)%%2==worker-1,,drop=FALSE]
    stopifnot(nrow(tasks)==1000,sum(tasks$replicate<=50)==200,sum(tasks$replicate>50)==800)
  }

  scratch = tempfile("extend_",tmpdir=file.path(root,"validation"))
  dir.create(scratch)
  scratch = normalizePath(scratch,winslash="/")
  old_plot = plot_study
  old_sim = get_sim
  plot_note = NULL
  assign("plot_study",function(summary,curves,calibration,cfg,fdir,note){
    plot_note <<- note
  },envir=.GlobalEnv)
  assign("get_sim",function(...) stop("Unexpected refit during extension check"),envir=.GlobalEnv)
  on.exit({
    assign("plot_study",old_plot,envir=.GlobalEnv)
    assign("get_sim",old_sim,envir=.GlobalEnv)
    allowed = paste0(normalizePath(file.path(root,"validation"),winslash="/"),"/")
    if (startsWith(scratch,allowed)) unlink(scratch,recursive=TRUE)
  })
  out = file.path(scratch,"results",cfg$run_id,"replicates")
  dir.create(out,recursive=TRUE)
  X = matrix(1,3,2,dimnames=list(NULL,c("SNP-1","SNP-2")))
  scored = score_fit(list(cs=list(1L),pip=c(.96,.01),converged=TRUE),"fSuSiE",X,1L,.9)
  row = cbind(data.frame(method="fSuSiE",success=TRUE),scored$stats,seconds=1,warnings="",error="")
  result = list(summary=row,PIP=scored$PIP,CS=list(fSuSiE=scored$CS))
  task = old_grid[1,,drop=FALSE]
  fixture = list(task=task,config=old_cfg,signature=old_sig,result=result,
    truth=list(n_snps=2L,pve_realized=.2),seconds=1)
  path = file.path(out,paste0(task$task_id,".rds"))
  atomic_save(fixture,path)
  before = tools::md5sum(path)
  stopifnot(run_one(scratch,cfg,grid[1,,drop=FALSE],sig)=="saved")
  duplicate = fixture; duplicate$signature = sig; duplicate$task = grid[1,,drop=FALSE]
  atomic_save(duplicate,path)
  stopifnot(identical(tools::md5sum(path),before),identical(readRDS(path),fixture))
  bad_task = grid[1,,drop=FALSE]; bad_task$seed = bad_task$seed+1
  stopifnot(inherits(tryCatch(run_one(scratch,cfg,bad_task,sig),error=identity),"error"))
  # The summary must combine old and newly extended checkpoints and use the
  # current target, not the 50-replicate setting stored in the first file.
  fixture$task = grid[which(grid$replicate==51)[1],,drop=FALSE]
  fixture$config = cfg; fixture$signature = sig
  atomic_save(fixture,file.path(out,paste0(fixture$task$task_id,".rds")))
  invisible(capture.output({ans_summary = summarize_study(scratch,cfg)}))
  stopifnot(ans_summary$n_completed==2,ans_summary$n_expected==2000,
    ans_summary$target_replicates==250,!ans_summary$complete,plot_note=="PARTIAL RUN",
    nrow(ans_summary$metrics)==2,ans_summary$summary$n_success==2)
  short = cfg; short$n_rep = 50L
  stopifnot(inherits(tryCatch(summarize_study(scratch,short),error=identity),"error"))
}
check_extension()
cat("PASS: legacy checkpoint reuse, unchanged seeds, 800 new datasets per worker,\n",
    "mixed-count summary at target 250, unchanged saved results, and rejection of incompatible changes.\n")
