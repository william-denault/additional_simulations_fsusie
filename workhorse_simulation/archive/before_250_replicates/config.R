# Change run_id when changing the simulation or fitting settings.
config = list(
  run_id="workhorse_50_P100_L12_T60_RV", n_rep=50L, seed=666L,
  L=1:2, f_type=c("step", "smooth"), noise_type=c("iid", "GP"),
  pve=.2, n=100L, n_pos=60L, n_snps=100L,
  gp_length_scale=12, gp_nugget=.15,
  L_fit=10L, coverage=.95, min_abs_corr=.5, pip_threshold=.9,
  pip_thresholds=sort(unique(c(seq(.01,.99,.02), .5,.8,.9,.95,.99,.999))),
  show_intervals=FALSE, fdr_plot_max=.25
)
