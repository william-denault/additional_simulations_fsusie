read_config = function(root){
  e = new.env(parent=baseenv())
  sys.source(file.path(root,"config.R"), envir=e)
  cfg = e$config
  stopifnot(cfg$n_rep>=1, cfg$n_rep==as.integer(cfg$n_rep),
    all(cfg$L>=1), all(cfg$L<=cfg$L_fit),
    grepl("^[A-Za-z0-9_-]+$",cfg$run_id))
  cfg
}

study_grid = function(cfg){
  g = expand.grid(L=cfg$L, f_type=cfg$f_type, noise_type=cfg$noise_type,
    pve=cfg$pve, replicate=seq_len(cfg$n_rep), stringsAsFactors=FALSE)
  g$task_id = with(g, paste0(f_type,"_",noise_type,"_L",L,"_pve",
    gsub("\\.","p",format(pve,trim=TRUE,scientific=FALSE)),"_rep",sprintf("%03d",replicate)))
  # Matching X, causal SNPs and block starts across shape/noise/PVE conditions.
  g$seed = cfg$seed+1000L*g$replicate+10L*g$L
  stopifnot(!anyDuplicated(g$task_id), all(is.finite(g$seed)))
  g
}

study_signature = function(root,cfg){
  fit_config = cfg[setdiff(names(cfg),c("show_intervals","fdr_plot_max","pip_thresholds"))]
  packages = c("fsusieR","mvsusieR","susieR","mashr")
  list(config=fit_config,
    packages=setNames(vapply(packages,function(p) as.character(utils::packageVersion(p)),""),packages),
    code=setNames(unname(tools::md5sum(file.path(root,"R",c("get_sim.R","run_methods.R")))),
                  c("get_sim.R","run_methods.R")))
}

atomic_save = function(x,path){
  tmp = tempfile("checkpoint_",tmpdir=dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp))
  saveRDS(x,tmp)
  # A duplicate worker may have finished the same deterministic task meanwhile.
  if (file.exists(path)){
    old = readRDS(path)
    if (!identical(old$signature,x$signature) || !identical(old$task,x$task))
      stop("A conflicting checkpoint already exists: ",path)
    return(invisible(FALSE))
  }
  if (!file.rename(tmp,path)){
    if (!file.exists(path) || !identical(readRDS(path)$signature,x$signature))
      stop("Could not save checkpoint: ",path)
  }
  invisible(TRUE)
}

run_one = function(root,cfg,task,signature){
  outdir = file.path(root,"results",cfg$run_id,"replicates")
  path = file.path(outdir,paste0(task$task_id,".rds"))
  if (file.exists(path)){
    old = readRDS(path)
    if (!identical(old$signature,signature))
      stop("Settings/packages/code differ from saved results. Choose a new run_id in config.R.")
    return("saved")
  }
  start = proc.time()[[3]]
  set.seed(task$seed)
  sim = get_sim(L=task$L,f_type=task$f_type,pve=task$pve,n=cfg$n,n_pos=cfg$n_pos,
    noise_type=task$noise_type,gp_length_scale=cfg$gp_length_scale,gp_nugget=cfg$gp_nugget)
  # The true V, E, F and noise_var cannot be accessed by the fitting functions.
  ans = run_methods(sim[c("X","Y","true_pos")],L_fit=cfg$L_fit,
    coverage=cfg$coverage,min_abs_corr=cfg$min_abs_corr,pip_threshold=cfg$pip_threshold)
  truth = list(true_pos=sim$true_pos,true_snps=colnames(sim$X)[sim$true_pos],
    n_snps=ncol(sim$X),starts=sim$starts,pve_realized=sim$pve_realized,
    signal_var=sim$signal_var,noise_var=sim$noise_var)
  atomic_save(list(task=task,config=cfg,signature=signature,result=ans,truth=truth,
    seconds=proc.time()[[3]]-start),path)
  "new"
}

run_study = function(root,cfg=read_config(root),worker=1L,n_workers=1L){
  stopifnot(worker>=1,worker<=n_workers)
  runroot = file.path(root,"results",cfg$run_id)
  for (d in c("replicates","logs")) dir.create(file.path(runroot,d),recursive=TRUE,showWarnings=FALSE)
  signature = study_signature(root,cfg)
  tasks = study_grid(cfg)
  tasks = tasks[(tasks$replicate-1L) %% n_workers == worker-1L,,drop=FALSE]
  log = file.path(runroot,"logs",paste0("worker_",worker,"_of_",n_workers,".txt"))
  writeLines(capture.output(sessionInfo()),file.path(runroot,"logs",paste0("session_worker_",worker,".txt")))
  say = function(message){
    line = paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),message)
    cat(line,"\n"); cat(line,"\n",file=log,append=TRUE)
    flush.console()
  }
  say(paste("Run",cfg$run_id,"worker",worker,"of",n_workers,"datasets",nrow(tasks)))
  counts = c(new=0L,saved=0L)
  for (i in seq_len(nrow(tasks))){
    status = run_one(root,cfg,tasks[i,,drop=FALSE],signature)
    counts[status] = counts[status]+1L
    say(sprintf("[%d/%d] %s: %s",i,nrow(tasks),tasks$task_id[i],status))
  }
  say(paste(names(counts),counts,collapse="; "))
  invisible(counts)
}

mean_se = function(x){
  x = x[is.finite(x)]
  n = length(x)
  avg = if (n) mean(x) else NA_real_
  se = if (n>1) sd(x)/sqrt(n) else NA_real_
  c(mean=avg,se=se,lower=pmax(0,avg-1.96*se),upper=pmin(1,avg+1.96*se))
}

group_rows = function(d,columns){
  split(seq_len(nrow(d)),interaction(d[columns],drop=TRUE,lex.order=TRUE))
}

condition_summary = function(metrics){
  keys = c("L","f_type","noise_type","pve","method")
  rows = lapply(group_rows(metrics,keys),function(ii){
    all = metrics[ii,,drop=FALSE]
    d = all[all$success,,drop=FALSE]
    row = all[1,keys,drop=FALSE]
    row$n_attempted = nrow(all)
    row$n_success = nrow(d)
    row$n_error = sum(!all$success)
    row$n_nonconverged = sum(d$converged %in% FALSE)
    row$n_unknown_convergence = sum(is.na(d$converged))
    row$n_with_cs = sum(d$n_cs>0)
    for (metric in c("coverage","cs_power","cs_available","pip_power","pip_fdp")){
      stat = mean_se(d[[metric]])
      for (name in names(stat)) row[[paste(metric,name,sep="_")]] = unname(stat[name])
    }
    total = sum(d$n_cs)
    row$pooled_coverage = if (total>0) sum(d$covered_cs)/total else NA_real_
    row$mean_cs_size = if (total>0) sum(d$total_cs_size)/total else NA_real_
    row$mean_n_cs = if (nrow(d)) mean(d$n_cs) else NA_real_
    row$mean_seconds = if (nrow(all)) mean(all$seconds) else NA_real_
    row
  })
  ans = do.call(rbind,rows)
  rownames(ans) = NULL
  ans
}

pip_curves = function(pips,thresholds){
  keys = c("task_id","replicate","L","f_type","noise_type","pve","method")
  rows = lapply(group_rows(pips,keys),function(ii){
    d = pips[ii,,drop=FALSE]
    if (!all(is.finite(d$PIP))) return(NULL)
    ans = d[rep(1,length(thresholds)),keys,drop=FALSE]
    ans$threshold = thresholds
    ans$power = vapply(thresholds,function(t) mean(d$PIP[d$causal]>=t),0.)
    ans$fdp = vapply(thresholds,function(t){
      selected = d$PIP>=t
      sum(selected & !d$causal)/max(1,sum(selected))
    },0.)
    ans
  })
  per_rep = do.call(rbind,rows)
  if (is.null(per_rep)) return(data.frame())
  keys = c("L","f_type","noise_type","pve","method","threshold")
  ans = lapply(group_rows(per_rep,keys),function(ii){
    d = per_rep[ii,,drop=FALSE]
    row = d[1,keys,drop=FALSE]
    row$n = nrow(d)
    row$power = mean(d$power)
    row$fdr = mean(d$fdp)
    row$power_se = if (nrow(d)>1) sd(d$power)/sqrt(nrow(d)) else NA_real_
    row$fdr_se = if (nrow(d)>1) sd(d$fdp)/sqrt(nrow(d)) else NA_real_
    row
  })
  do.call(rbind,ans)
}

pip_calibration = function(pips){
  pips = pips[is.finite(pips$PIP),,drop=FALSE]
  if (!nrow(pips)) return(data.frame())
  keys = c("f_type","noise_type","pve","method")
  rows = lapply(group_rows(pips,keys),function(ii){
    d = pips[ii,,drop=FALSE]
    # Pool L, but cluster the uncertainty by replicate (paired conditions).
    reps = sort(unique(d$replicate))
    bin = pmin(10L,as.integer(d$PIP*10)+1L)
    bins = lapply(1:10,function(b){
      use = bin==b
      if (!any(use)) return(NULL)
      ns = vapply(reps,function(r) sum(use & d$replicate==r),0.)
      ys = vapply(reps,function(r) sum(use & d$replicate==r & d$causal),0.)
      frequency = sum(ys)/sum(ns)
      se = if (length(reps)>1) sd(ys-frequency*ns)/sqrt(length(reps))/mean(ns) else NA_real_
      cbind(d[1,keys,drop=FALSE],bin=b,n_snps=sum(use),mean_pip=mean(d$PIP[use]),
        observed=frequency,se=se,lower=max(0,frequency-1.96*se),upper=min(1,frequency+1.96*se))
    })
    do.call(rbind,bins)
  })
  do.call(rbind,rows)
}

summarize_study = function(root,cfg=read_config(root)){
  paths = list.files(file.path(root,"results",cfg$run_id,"replicates"),"\\.rds$",full.names=TRUE)
  if (!length(paths)) stop("No saved simulations for ",cfg$run_id)
  records = lapply(paths,readRDS)
  signature = records[[1]]$signature
  if (!all(vapply(records,function(x) identical(x$signature,signature),logical(1))))
    stop("The saved run mixes settings, code or package versions.")
  plot_config = cfg[c("show_intervals","fdr_plot_max","pip_thresholds")]
  cfg = records[[1]]$config
  cfg[names(plot_config)] = plot_config
  expected = study_grid(cfg)
  ids = vapply(records,function(x) x$task$task_id,"")
  stopifnot(!anyDuplicated(ids),all(ids %in% expected$task_id))
  metrics = do.call(rbind,lapply(records,function(x){
    d = x$result$summary
    cbind(x$task[rep(1,nrow(d)),,drop=FALSE],d,
      pve_realized=x$truth$pve_realized,n_snps=x$truth$n_snps)
  }))
  pips = do.call(rbind,lapply(records,function(x){
    d = x$result$PIP
    cbind(x$task[rep(1,nrow(d)),,drop=FALSE],d)
  }))
  rownames(metrics) = rownames(pips) = NULL
  summary = condition_summary(metrics)
  curves = pip_curves(pips,cfg$pip_thresholds)
  calibration = pip_calibration(pips)
  tdir = file.path(root,"tables",cfg$run_id)
  fdir = file.path(root,"figures",cfg$run_id)
  dir.create(tdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(fdir,recursive=TRUE,showWarnings=FALSE)
  outputs = list(metrics_per_simulation=metrics,PIPs_per_simulation=pips,
    summary_by_condition=summary,PIP_power_FDR=curves,PIP_calibration=calibration)
  converged = metrics[metrics$success & metrics$converged %in% TRUE,,drop=FALSE]
  if (nrow(converged)) outputs$summary_converged_only = condition_summary(converged)
  for (name in names(outputs)) utils::write.csv(outputs[[name]],file.path(tdir,paste0(name,".csv")),row.names=FALSE)
  complete = length(records)==nrow(expected)
  note = if (grepl("pilot",cfg$run_id)) {
    if (complete) "PILOT" else "PARTIAL PILOT"
  } else if (!complete) "PARTIAL RUN" else ""
  plot_study(summary,curves,calibration,cfg,fdir,note)
  minutes = sum(vapply(records,function(x) x$seconds,0.))/60
  report = c(paste("Run:",cfg$run_id),
    sprintf("Completed datasets: %d / %d. Replicates per condition: %d.",length(records),nrow(expected),cfg$n_rep),
    sprintf("Fit errors: %d. Nonconverged fits: %d. Unknown convergence: %d.",
      sum(!metrics$success),sum(metrics$success & metrics$converged %in% FALSE),
      sum(metrics$success & is.na(metrics$converged))),
    sprintf("Dataset runtime sum: %.1f min; projected full sequential run: %.1f min; ideal two workers: %.1f min.",
      minutes,minutes*nrow(expected)/length(records),minutes*nrow(expected)/length(records)/2),
    "All methods estimate noise from observed data; no true residual covariance is passed to a fit.",
    paste0("fSuSiE uses standardize_Y=TRUE. All methods use L_fit=",cfg$L_fit," (true L is not supplied)."),
    "Coverage: mean per-simulation fraction of CSs containing a causal SNP; simulations with no CS are excluded.",
    "CS power: mean fraction of causal SNPs in the union of reported CSs; no CS gives power zero.",
    "PIP FDR: mean of FP/max(1,discoveries) across simulations, separately from CS coverage/power.",
    "Failed fits are excluded from accuracy metrics; nonconverged fits remain included and flagged.",
    "Mean +/- 1.96 Monte Carlo SE is exported; intervals are approximate and coarse near boundaries.",
    "Check CS availability/size, convergence and failures alongside coverage.",
    paste("Packages:",paste(names(signature$packages),signature$packages,collapse="; ")),
    paste("Tables:",normalizePath(tdir,winslash="/")),paste("Figures:",normalizePath(fdir,winslash="/")))
  writeLines(report,file.path(tdir,"run_report.txt"))
  cat(paste(report,collapse="\n"),"\n")
  invisible(list(summary=summary,metrics=metrics,curves=curves,calibration=calibration))
}
