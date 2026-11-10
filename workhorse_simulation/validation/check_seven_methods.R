# Run from the project folder. Two execution checks, not a performance study.
root = normalizePath(getwd(),winslash="/")
source(file.path(root,"R","setup.R"))
setup(root)
cfg = read_config(root)
cfg$n_rep = 1L
cfg$run_id = paste0(cfg$run_id,"_validation")
for (folder in c("replicates","logs"))
  dir.create(file.path(root,"results",cfg$run_id,folder),recursive=TRUE,showWarnings=FALSE)
grid = study_grid(cfg)
tasks = grid[(grid$f_type=="step" & grid$noise_type=="GP" & grid$L==2) |
               (grid$f_type=="smooth" & grid$noise_type=="iid" & grid$L==1),,drop=FALSE]
signature = study_signature(root,cfg)
for (i in seq_len(nrow(tasks))){
  task = tasks[i,,drop=FALSE]
  cat("Checking",task$task_id,"with seven methods.\n"); flush.console()
  run_one(root,cfg,task,signature)
  path = file.path(root,"results",cfg$run_id,"replicates",paste0(task$task_id,".rds"))
  record = readRDS(path)
  ans = record$result
  print(ans$summary[,c("method","success","converged","residual_change_from_initial","seconds","error")])
  fixed = ans$summary$residual_policy=="fixed"
  updated = ans$summary$residual_policy=="estimated"
  stopifnot(nrow(ans$summary)==7,all(ans$summary$success),record$truth$n_snps<=100,
    nrow(ans$PIP)==7*record$truth$n_snps,all(is.finite(ans$PIP$PIP)),
    sum(fixed)==3,sum(updated)==3,
    all(ans$summary$residual_change_from_initial[fixed]<1e-8),
    all(is.finite(ans$summary$residual_change_from_initial[updated])),
    any(ans$summary$residual_change_from_initial[updated]>1e-8))
  for (method in names(ans$CS)){
    sets = ans$CS[[method]]
    row = ans$summary[ans$summary$method==method,,drop=FALSE]
    truth = record$truth$true_pos
    coverage = if (length(sets)) mean(vapply(sets,function(s) any(s %in% truth),TRUE)) else NA_real_
    stopifnot(isTRUE(all.equal(row$coverage,coverage)),
      row$cs_power==length(intersect(unlist(sets),truth))/length(truth))
  }
  stopifnot(run_one(root,cfg,task,signature)=="saved")
}
summarize_study(root,cfg)
cat("PASS: seven fits per dataset, fixed/updated residual covariance, metrics, restart and plotting.\n")
