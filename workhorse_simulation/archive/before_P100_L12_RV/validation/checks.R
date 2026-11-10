# Run from the project folder: Rscript --vanilla validation/checks.R
root = normalizePath(getwd(),winslash="/")
source(file.path(root,"R","setup.R"))
setup(root)
cfg = read_config(root)

set.seed(24)
a = get_sim(L=3,f_type="step",pve=.2,n_pos=60,noise_type="iid")
set.seed(24)
b = get_sim(L=3,f_type="step",pve=.2,n_pos=60,noise_type="GP")
stopifnot(identical(a$X,b$X),identical(a$F,b$F),
  identical(dim(b$Y),c(100L,60L)),
  setequal(which(rowSums(abs(b$F))>0),b$true_pos),
  all(rowSums(b$F[b$true_pos,,drop=FALSE]!=0)==40),
  isTRUE(all.equal(b$Y,b$X%*%b$F+b$E)),
  isTRUE(all.equal(b$signal_var/(b$signal_var+mean(diag(b$V))),.2)),
  isTRUE(all.equal(diag(b$V),rep(b$noise_var,60),check.attributes=FALSE)),
  abs(b$V[1,13]/b$noise_var-.85*exp(-.5))<1e-12)
set.seed(24)
c = get_sim(L=3,f_type="smooth",pve=.2,n_pos=60,noise_type="GP")
stopifnot(identical(a$X,c$X),identical(a$true_pos,c$true_pos),
  all((rowSums(abs(c$F))>0) == (seq_len(nrow(c$F)) %in% c$true_pos)))

X = matrix(1,5,4,dimnames=list(NULL,paste0("SNP-",1:4)))
fit = list(cs=list(c(1),integer(),c(4),c(1)),pip=c(.96,.01,.92,.001),converged=TRUE)
s = score_fit(fit,"fSuSiE",X,c(1,3),.9)
stopifnot(s$stats$n_cs==2,s$stats$coverage==.5,s$stats$cs_power==.5,
          s$stats$pip_power==1,s$stats$pip_fdp==0)
fit$cs = list()
s0 = score_fit(fit,"fSuSiE",X,c(1,3),.9)
stopifnot(is.na(s0$stats$coverage),s0$stats$cs_power==0,s0$stats$n_cs==0)
fit$cs = list(1:2)
s1 = score_fit(fit,"fSuSiE",X,c(1,3),.9)
d = rbind(s$stats,s1$stats,s0$stats)
d$L = 2; d$f_type = "step"; d$noise_type = "iid"; d$pve = .2
d$method = "fSuSiE"; d$success = TRUE; d$seconds = 1
summ = condition_summary(d)
stopifnot(summ$coverage_mean==.75,summ$pooled_coverage==2/3,
          summ$cs_power_mean==1/3,summ$n_with_cs==2)
grid = study_grid(cfg)
odd = grid$task_id[(grid$replicate-1)%%2==0]
even = grid$task_id[(grid$replicate-1)%%2==1]
stopifnot(nrow(grid)==600,length(intersect(odd,even))==0,
          setequal(c(odd,even),grid$task_id))
check_restart = function(){
  scratch = tempfile("restart_",tmpdir=file.path(root,"validation"))
  dir.create(scratch)
  scratch = normalizePath(scratch,winslash="/")
  on.exit({
    allowed = paste0(normalizePath(file.path(root,"validation"),winslash="/"),"/")
    if (startsWith(scratch,allowed)) unlink(scratch,recursive=TRUE)
  })
  path = file.path(scratch,"results",cfg$run_id,"replicates",paste0(grid$task_id[1],".rds"))
  dir.create(dirname(path),recursive=TRUE)
  task = grid[1,,drop=FALSE]
  signature = study_signature(root,cfg)
  fixture = list(task=task,signature=signature)
  atomic_save(fixture,path)
  before = file.info(path)$mtime
  stopifnot(run_one(scratch,cfg,task,signature)=="saved",
    identical(before,file.info(path)$mtime))
  bad = signature; bad$config$seed = 999L
  stopifnot(inherits(tryCatch(run_one(scratch,cfg,task,bad),error=identity),"error"))
  atomic_save(fixture,path)
  stopifnot(identical(readRDS(path),fixture))
}
check_restart()
warning_check = capture_fit(function(){
  cat("WARNING: ELBO decreased in a test fixture\n")
  warning("test R warning"); 123
})
stopifnot(warning_check$success,grepl("ELBO decreased",warning_check$warnings),
          grepl("test R warning",warning_check$warnings))
cat("PASS: generation, GP covariance, causal assignment, zero-CS handling,\n",
    "simulation-wise coverage, CS power, worker partition, restart and warning capture.\n")
