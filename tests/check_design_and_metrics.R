# Run from the project directory: source("tests/check_design_and_metrics.R")
source("R/study.R")
source("R/summarize.R")
cfg <- load_config()
validate_config(cfg)
f_method <- method_grid(cfg, "gp")
f_method <- f_method[f_method$engine == "fSuSiE", ]
stopifnot(nrow(f_method) == 1L,
  f_method$method == if (isTRUE(cfg$fsusie_standardize_Y)) "fSuSiE_standardized" else "fSuSiE_default")
tasks <- make_tasks(cfg)
stopifnot(nrow(tasks) == 12L * cfg$n_rep, !anyDuplicated(tasks$task_id))
# The two background workers partition all tasks, without overlap.
a <- tasks$task_id[(tasks$replicate - 1L) %% 2L == 0L]
b <- tasks$task_id[(tasks$replicate - 1L) %% 2L == 1L]
stopifnot(!length(intersect(a, b)), setequal(c(a, b), tasks$task_id))

step <- effect_profiles(cfg, "step")
smooth <- effect_profiles(cfg, "smooth")
stopifnot(isTRUE(all.equal(rowMeans(step^2), rep(cfg$effect_mean_square, 3))),
          isTRUE(all.equal(rowMeans(smooth^2), rep(cfg$effect_mean_square, 3))),
          identical(step != 0, smooth != 0))
v <- noise_covariance(cfg, "gp")
stopifnot(isTRUE(all.equal(diag(v), rep(cfg$noise_sd^2, cfg$positions))),
          isTRUE(all.equal(crossprod(chol(v)), v)))
d1 <- simulate_dataset(cfg, tasks[1, ])
d2 <- simulate_dataset(cfg, tasks[2, ])
stopifnot(identical(d1$X, d2$X), identical(d1$V, d2$V),
          identical(d1$permutation, d2$permutation),
          isTRUE(all.equal(d1$Y - d1$X[, d1$causal, drop = FALSE] %*% d1$B,
                           d2$Y - d2$X[, d2$causal, drop = FALSE] %*% d2$B)),
          identical(d1, simulate_dataset(cfg, tasks[1, ])))

z <- score_pips(c(0.95, 0.91, 0.2, 0.05), c(1L, 3L), c(0.9, 0.99))
stopifnot(z$tp[1] == 1, z$fp[1] == 1, z$power[1] == 0.5, z$fdp[1] == 0.5,
          z$discoveries[2] == 0, z$power[2] == 0, z$fdp[2] == 0,
          identical(credible_set(c(0.6, 0.3, 0.1), 0.9), 1:2))
z <- score_sets(list(c(1L, 2L), 4L), c(1L, 3L))
stopifnot(z['n_cs'] == 2, z['covered_cs'] == 1, z['cs_power'] == 0.5,
          is.na(ratio_summary(c(0, 0), c(0, 0), 1:2)[1]),
          ratio_summary(c(1, 2), c(1, 3), 1:2)[1] == 0.75)
x <- data.frame(method = 'test', threshold = 0.9, replicate = 1:2,
                 power = c(1, 1), fdp = c(0, 0.75), tp = c(1, 1), fp = c(0, 3),
                 discoveries = c(1, 4))
z <- aggregate_pips(x, "method")
stopifnot(z$fdr == 0.375, z$pooled_fdr == 0.6)
# Coverage excludes zero-CS fits; CS power must keep those as zero power.
x <- data.frame(method = "test", nominal = 0.95, replicate = 1:3,
  n_cs = c(0, 1, 3), covered_cs = c(0, 1, 1), total_size = c(0, 2, 6),
  cs_power = c(0, 0.5, 1))
z <- aggregate_coverage(x, "method")
stopifnot(z$n_success == 3, z$n_with_cs == 2, z$cs_coverage == 0.5,
  z$cs_power == 0.5, z$total_cs == 4,
  isTRUE(all.equal(z$coverage_mcse, unname(ratio_summary(c(1, 1), c(1, 3), 2:3)[2]))))
z <- aggregate_coverage(x[1, ], "method")
stopifnot(is.na(z$cs_coverage), is.na(z$coverage_mcse), z$n_with_cs == 0,
  z$cs_power == 0)
z <- score_sets(list(c(1L, 2L), c(1L, 4L)), c(1L, 3L))
stopifnot(z["cs_power"] == 0.5) # A causal SNP in two CSs is counted only once.
cat("PASS: paired design, matched effects, GP covariance, worker partition, CS coverage, and power/FDR definitions.\n")
