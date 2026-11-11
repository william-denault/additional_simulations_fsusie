# Edit this file before starting either worker. Use a new run_id after edits.
study_config <- list(
  run_id = "main_100_T60_gpL12_stdY",
  seed = 20261110L,
  n_rep = 100L,
  effect_types = c("step", "smooth"),
  noise_types = c("iid", "gp"),
  true_L = 1:3,
  n = 100L,
  p = 12L,                   # Original step vignette size; increase to 100 later.
  positions = 60L,           # Divisible by 4; fSuSiE handles its own wavelet remapping.
  fit_L = 5L,                 # SAME upper limit for all fits; true L varies.
  fsusie_standardize_Y = TRUE, # Scale phenotype positions AND wavelet coefficients.
  effect_mean_square = 0.10,  # Per causal SNP, averaged over ALL positions.
  noise_sd = 1,
  gp_length_fraction = 0.20,  # Length scale = 0.20 * 60 = 12 positions.
  gp_nugget = 0.15,          # V = sd^2 * ((1-nugget)*RBF + nugget*I).
  proxy_copy_probability = 0.98,
  priors = c("diagonal", "canonical", "mash"),
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
  include_shuffled_fsusie = FALSE, # Only default fSuSiE and the mvSuSiE comparisons.
  save_full_fits = FALSE
)
