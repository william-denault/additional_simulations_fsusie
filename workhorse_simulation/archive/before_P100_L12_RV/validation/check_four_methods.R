# Run from the project folder. One actual GP dataset; not a performance study.
root = normalizePath(getwd(),winslash="/")
source(file.path(root,"R","setup.R"))
setup(root)
set.seed(666)
sim = get_sim(L=2,f_type="step",pve=.2,n=100,n_pos=60,noise_type="GP")
observed = sim[c("X","Y","true_pos")]
stopifnot(is.null(observed$V),is.null(observed$noise_var),is.null(observed$E))
cat("Checking all four methods on one GP dataset, without supplying V.\n")
ans = run_methods(observed)
print(ans$summary[,c("method","success","n_cs","coverage","cs_power","converged","seconds","error")])
stopifnot(nrow(ans$summary)==4,all(ans$summary$success),
          nrow(ans$PIP)==4*ncol(sim$X),all(is.finite(ans$PIP$PIP)))
for (method in names(ans$CS)){
  sets = ans$CS[[method]]
  row = ans$summary[ans$summary$method==method,,drop=FALSE]
  coverage = if (length(sets)) mean(vapply(sets,function(s) any(s %in% sim$true_pos),TRUE)) else NA_real_
  stopifnot(isTRUE(all.equal(row$coverage,coverage)),
    row$cs_power==length(intersect(unlist(sets),sim$true_pos))/length(sim$true_pos))
}
saveRDS(ans,file.path(root,"validation","GP_check.rds"))
cat("PASS: four real GP fits, per-simulation coverage, causal SNP capture and PIPs.\n")
