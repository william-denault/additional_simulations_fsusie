# Edit this file before starting either worker. Use a new run_id after edits.
study_config <- list(
  run_id = "main_50",
  seed = 20261110L,
  n_rep = 50L,
  effect_types = c("step", "smooth"),
  noise_types = c("iid", "gp"),
  true_L = 1:3,
  n = 100L,
  p = 12L,                   # Original step vignette size; increase to 100 later.
  positions = 32L,             # Power of two; divisible by 4. 64/128 cost more.
  fit_L = 5L,                 # SAME upper limit for all fits; true L varies.
  effect_mean_square = 0.10,  # Per causal SNP, averaged over ALL positions.
  noise_sd = 1,
  gp_length_fraction = 0.10,  # Length scale = this fraction * positions.
  gp_nugget = 0.15,          # V = sd^2 * ((1-nugget)*RBF + nugget*I).
  proxy_copy_probability = 0.98,
  priors = c("diagonal", "canonical", "mash", "oracle_shape"),
  residuals = c("true_fixed", "estimated_empirical_init", "estimated_true_init"),
  # Add "empirical_fixed" above to reproduce the fourth initialization policy.
  estimate_prior_variance = TRUE,
  # Freeze each supplied mixture across covariance policies for a controlled
  # comparison. Effect-specific overall scales are still estimated above.
  estimate_prior_mixture_weights = FALSE,
  mv_max_iter = 100L,
  mv_tol = 1e-4,
  coverage = 0.95,
  min_abs_corr = 0.5,
  coverage_levels = c(0.50, 0.70, 0.80, 0.90, 0.95, 0.99),
  pip_thresholds = sort(unique(c(seq(0.01, 0.99, by = 0.02),
                                 0.50, 0.80, 0.90, 0.95, 0.99, 0.999))),
  table_pip_threshold = 0.90,
  include_shuffled_fsusie = TRUE,  # IID ONLY: controls the residual distribution.
  save_full_fits = FALSE
)
