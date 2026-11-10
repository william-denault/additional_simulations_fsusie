# Simulation, fitting, scoring, and resumable execution. No package installation.

load_config <- function(root = getwd()) {
  e <- new.env(parent = baseenv())
  sys.source(file.path(root, "config", "study_config.R"), envir = e)
  e$study_config
}

check_dependencies <- function() {
  # Rscript on this laptop does not automatically include the existing user library.
  if (.Platform$OS.type == "windows") {
    v <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
    lib <- c(file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", v),
             file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", v))
    .libPaths(unique(c(.libPaths(), lib[dir.exists(lib)])))
  }
  pkgs <- c("fsusieR", "mvsusieR", "susieR", "mashr", "ggplot2")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the missing packages first: ", paste(missing, collapse = ", "))
  versions <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
  print(versions)
  invisible(versions)
}

validate_config <- function(cfg) {
  positive_int <- c("n_rep", "n", "p", "positions", "fit_L", "mv_max_iter", "seed")
  stopifnot(all(vapply(cfg[positive_int], function(x)
    length(x) == 1L && is.finite(x) && x >= 1 && x == as.integer(x), logical(1))))
  stopifnot(cfg$n >= 20L, cfg$p >= 12L, cfg$positions >= 16L,
            cfg$positions %% 4L == 0L, # Equal quarter-block effects; fSuSiE remaps internally.
            all(cfg$true_L %in% 1:3), cfg$fit_L >= max(cfg$true_L),
            all(cfg$effect_types %in% c("step", "smooth")),
            all(cfg$noise_types %in% c("iid", "gp")),
            all(cfg$priors %in% c("diagonal", "canonical", "mash", "oracle_shape")),
            all(cfg$residuals %in% c("true_fixed", "estimated_empirical_init",
                                     "estimated_true_init", "empirical_fixed")),
            cfg$effect_mean_square > 0, cfg$noise_sd > 0,
            cfg$gp_length_fraction > 0, cfg$gp_nugget > 0, cfg$gp_nugget <= 1,
            cfg$proxy_copy_probability >= 0, cfg$proxy_copy_probability <= 1,
            all(cfg$coverage_levels > 0 & cfg$coverage_levels < 1),
            all(cfg$pip_thresholds > 0 & cfg$pip_thresholds <= 1),
            cfg$table_pip_threshold %in% cfg$pip_thresholds,
            cfg$coverage %in% cfg$coverage_levels,
            grepl("^[A-Za-z0-9_-]+$", cfg$run_id))
  invisible(TRUE)
}

make_tasks <- function(cfg) {
  x <- expand.grid(effect_type = cfg$effect_types, noise_type = cfg$noise_types,
                   true_L = cfg$true_L, replicate = seq_len(cfg$n_rep),
                   stringsAsFactors = FALSE)
  x$task_id <- sprintf("%s_%s_L%d_rep%03d", x$effect_type, x$noise_type,
                       x$true_L, x$replicate)
  x
}

effect_profiles <- function(cfg, shape) {
  t <- cfg$positions
  q <- t %/% 4L
  # Same supports and total squared effect for both shapes. The cosine bump is
  # the positive half-wave construction in the supplied smooth-effect vignette,
  # rescaled to the blocks used in the step vignette.
  block <- if (shape == "step") rep(1, q) else sin(pi * ((seq_len(q) - 0.5) / q))
  block <- block * sqrt(t * cfg$effect_mean_square / sum(block^2))
  b <- matrix(0, 3L, t)
  b[1, (q + 1L):(2L * q)] <- block
  b[2, (3L * q + 1L):t] <- block
  b[3, (q + 1L):(2L * q)] <- -block
  b
}

noise_covariance <- function(cfg, noise) {
  t <- cfg$positions
  if (noise == "iid") return(cfg$noise_sd^2 * diag(t))
  d <- outer(seq_len(t), seq_len(t), "-")
  kernel <- exp(-0.5 * (d / (cfg$gp_length_fraction * t))^2)
  cfg$noise_sd^2 * ((1 - cfg$gp_nugget) * kernel + cfg$gp_nugget * diag(t))
}

simulate_dataset <- function(cfg, task) {
  # Independent streams: paired X and standard-normal innovations across shape,
  # noise, and L for each replicate. No dependence on worker count or job order.
  set.seed(cfg$seed + 100L * task$replicate)
  maf <- runif(cfg$p, 0.10, 0.45)
  x <- matrix(rbinom(cfg$n * cfg$p, 2, rep(maf, each = cfg$n)), cfg$n, cfg$p)
  for (j in which(apply(x, 2, stats::var) == 0)) {
    repeat {
      x[, j] <- rbinom(cfg$n, 2, maf[j])
      if (stats::var(x[, j]) > 0) break
    }
  }
  # An actual allele-count proxy, rather than adding Gaussian noise to genotypes.
  x[, 4] <- x[, 3]
  replace <- runif(cfg$n) > cfg$proxy_copy_probability
  x[replace, 4] <- rbinom(sum(replace), 2, maf[3])
  if (stats::var(x[, 4]) == 0) x[, 4] <- x[, 3]
  x <- scale(x)
  colnames(x) <- paste0("SNP-", seq_len(cfg$p))
  causal <- c(1L, 3L, 9L)[seq_len(task$true_L)]
  b <- effect_profiles(cfg, task$effect_type)[seq_len(task$true_L), , drop = FALSE]
  v <- noise_covariance(cfg, task$noise_type)
  set.seed(cfg$seed + 100L * task$replicate + 1L)
  z <- matrix(rnorm(cfg$n * cfg$positions), cfg$n, cfg$positions)
  # chol(V) is upper triangular; Z %*% chol(V) has row covariance V.
  e <- z %*% chol(v)
  y <- x[, causal, drop = FALSE] %*% b + e
  colnames(y) <- paste0("CpG-", seq_len(cfg$positions))
  set.seed(cfg$seed + 100L * task$replicate + 2L)
  permutation <- sample.int(cfg$positions)
  list(X = x, Y = y, B = b, causal = causal, V = v, permutation = permutation,
       proxy_correlation = cor(x[, 3], x[, 4]))
}

method_grid <- function(cfg, noise) {
  m <- expand.grid(prior = cfg$priors, residual = cfg$residuals, stringsAsFactors = FALSE)
  m$method <- paste("mvSuSiE", m$prior, m$residual, sep = "__")
  m$engine <- "mvSuSiE"
  f <- data.frame(prior = "functional_default", residual = "fsusie_default",
                   method = "fSuSiE_default", engine = "fSuSiE")
  m <- rbind(m, f)
  if (cfg$include_shuffled_fsusie && noise == "iid") {
    f$method <- "fSuSiE_shuffled"
    m <- rbind(m, f)
  }
  m
}

capture_run <- function(fun) {
  warnings <- character()
  started <- proc.time()[[3]]
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(ok = !inherits(value, "error"), value = value,
       seconds = proc.time()[[3]] - started, warnings = paste(warnings, collapse = " | "))
}

make_prior <- function(kind, dat) {
  t <- ncol(dat$Y)
  if (kind == "canonical") return(mvsusieR::create_mixture_prior(R = t))
  if (kind == "diagonal") {
    return(mvsusieR::create_mixture_prior(mixture_prior =
      list(matrices = list(0.2 * diag(apply(dat$Y, 2, var))), weights = 1)))
  }
  if (kind == "oracle_shape") {
    # Oracle DIRECTION library, not oracle SNP identities. Symmetric normal
    # priors make opposite signs the same pattern; remove duplicate matrices.
    u <- lapply(seq_len(nrow(dat$B)), function(i) tcrossprod(dat$B[i, ]))
    key <- vapply(u, function(a) paste(signif(as.vector(a), 10), collapse = ","), character(1))
    u <- u[!duplicated(key)]
    return(mvsusieR::create_mixture_prior(mixture_prior =
      list(matrices = u, weights = rep(1 / length(u), length(u)))))
  }
  z <- susieR::calc_z(dat$X, dat$Y, center = TRUE, scale = TRUE)
  d <- mashr::mash_set_data(Bhat = z, Shat = matrix(1, nrow(z), ncol(z)))
  # Same simplified canonical-pattern recipe as the supplied vignette. It
  # learns weights/scales, not new spatial covariance patterns. Cache once per Y.
  invisible(utils::capture.output(
    g <- mashr::mash(d, Ulist = mashr::cov_canonical(d), outputlevel = 0)))
  mvsusieR::create_mixture_prior(fitted_g = g$fitted_g)
}

fit_method <- function(cfg, dat, method, prior = NULL) {
  if (method$engine == "fSuSiE") {
    y <- if (method$method == "fSuSiE_shuffled") dat$Y[, dat$permutation] else dat$Y
    # Preserve default fSuSiE prior, covariance estimation, postprocessing,
    # filtering, and optimization. Only set the shared L and output criteria.
    return(fsusieR::susiF(Y = y, X = dat$X, L = cfg$fit_L, verbose = FALSE,
                         cov_lev = cfg$coverage, min_purity = cfg$min_abs_corr))
  }
  args <- list(X = dat$X, Y = dat$Y, L = cfg$fit_L, prior_variance = prior,
    # X already has variance one. FALSE prevents mvSuSiE's extra rescaling of
    # the supplied effect covariance by observed phenotype SDs (oracle shape).
    standardize = FALSE, intercept = TRUE,
    estimate_prior_variance = cfg$estimate_prior_variance,
    estimate_prior_mixture_weights = cfg$estimate_prior_mixture_weights,
    max_iter = cfg$mv_max_iter, tol = cfg$mv_tol, verbose = FALSE, n_thread = 1L,
    coverage = cfg$coverage, min_abs_corr = cfg$min_abs_corr)
  policy <- method$residual
  args$estimate_residual_variance <- policy %in% c("estimated_empirical_init", "estimated_true_init")
  if (policy %in% c("true_fixed", "estimated_true_init")) args$residual_variance <- dat$V
  # For the empirical policies, NULL invokes the package's cov(Y) initialization.
  do.call(mvsusieR::mvsusie, args)
}

credible_set <- function(alpha, level) {
  stopifnot(all(is.finite(alpha)), all(alpha >= 0), abs(sum(alpha) - 1) < 1e-6)
  o <- order(alpha, decreasing = TRUE)
  k <- which(cumsum(alpha[o]) >= level - 1e-12)[1]
  if (is.na(k)) k <- length(o)
  sort(o[seq_len(k)])
}

posterior_summary <- function(fit, engine, p) {
  pip <- as.numeric(fit$pip)
  stopifnot(length(pip) == p, all(is.finite(pip)), all(pip >= -1e-10 & pip <= 1 + 1e-10))
  if (engine == "fSuSiE") {
    keep <- which(lengths(fit$cs) > 0)
    cs <- fit$cs[keep]
    a <- if (length(keep)) do.call(rbind, fit$alpha[keep]) else matrix(numeric(), 0, p)
  } else {
    cs <- fit$sets$cs
    if (length(cs)) {
      keep <- fit$sets$cs_index
      if (is.null(keep)) keep <- as.integer(sub("^L", "", names(cs)))
      stopifnot(length(keep) == length(cs), !anyNA(keep))
      a <- fit$alpha[keep, , drop = FALSE]
    } else a <- matrix(numeric(), 0, p)
  }
  # Drop identical reported sets consistently for scoring both methods.
  if (length(cs)) {
    keep <- !duplicated(vapply(cs, function(s) paste(sort(s), collapse = ","), character(1)))
    cs <- cs[keep]
    a <- a[keep, , drop = FALSE]
    stopifnot(ncol(a) == p, all(abs(rowSums(a) - 1) < 1e-6))
  }
  list(pip = pmin(1, pmax(0, pip)), cs = cs, selected_alpha = a,
       converged = if (is.null(fit$converged)) NA else isTRUE(fit$converged),
       niter = if (is.null(fit$niter)) NA_integer_ else as.integer(fit$niter))
}

score_pips <- function(pip, causal, thresholds) {
  do.call(rbind, lapply(thresholds, function(t) {
    selected <- which(pip >= t)
    tp <- sum(selected %in% causal)
    fp <- length(selected) - tp
    data.frame(threshold = t, tp = tp, fp = fp, discoveries = length(selected),
               power = tp / length(causal), fdp = fp / max(1L, length(selected)))
  }))
}

score_sets <- function(sets, causal) {
  hit <- vapply(sets, function(s) any(s %in% causal), logical(1))
  c(n_cs = length(sets), covered_cs = sum(hit),
    total_size = sum(lengths(sets)),
    recovered_causal = sum(causal %in% unique(unlist(sets))),
    cs_power = sum(causal %in% unique(unlist(sets))) / length(causal))
}

score_coverage <- function(post, causal, cfg) {
  # Freeze components reported at 95%; recompute membership only. This avoids
  # selecting different components at each nominal level. Report conditional
  # coverage and report CS availability separately (no CS => NA coverage).
  do.call(rbind, lapply(cfg$coverage_levels, function(level) {
    sets <- if (abs(level - cfg$coverage) < 1e-10) post$cs else
      lapply(seq_len(nrow(post$selected_alpha)), function(i)
        credible_set(post$selected_alpha[i, ], level))
    as.data.frame(as.list(c(nominal = level, score_sets(sets, causal))))
  }))
}

run_dataset <- function(cfg, task, signature) {
  began <- proc.time()[[3]]
  dat <- simulate_dataset(cfg, task)
  priors <- setNames(lapply(cfg$priors, function(k) capture_run(function() make_prior(k, dat))), cfg$priors)
  methods <- method_grid(cfg, task$noise_type)
  rows <- vector("list", nrow(methods))
  outputs <- vector("list", nrow(methods))
  for (i in seq_len(nrow(methods))) {
    m <- methods[i, ]
    prior <- if (m$engine == "mvSuSiE") priors[[m$prior]] else NULL
    if (!is.null(prior) && !prior$ok) {
      res <- prior
    } else {
      res <- capture_run(function() {
        fit <- fit_method(cfg, dat, m, if (!is.null(prior)) prior$value else NULL)
        post <- posterior_summary(fit, m$engine, cfg$p)
        list(post = post,
             pips = score_pips(post$pip, dat$causal, cfg$pip_thresholds),
             coverage = score_coverage(post, dat$causal, cfg),
             fit = if (cfg$save_full_fits) fit else NULL)
      })
    }
    row <- cbind(task, m, success = res$ok,
                 converged = if (res$ok) res$value$post$converged else NA,
                 niter = if (res$ok) res$value$post$niter else NA_integer_,
                 fit_seconds = res$seconds,
                 prior_seconds = if (is.null(prior)) 0 else prior$seconds,
                 warnings = paste(c(if (!is.null(prior)) prior$warnings, res$warnings), collapse = " | "),
                 error = if (res$ok) "" else conditionMessage(res$value))
    rows[[i]] <- row
    if (res$ok) outputs[[i]] <- res$value
  }
  names(outputs) <- methods$method
  list(signature = signature, task = task, status = do.call(rbind, rows), outputs = outputs,
       causal = dat$causal, proxy_correlation = dat$proxy_correlation,
       elapsed_seconds = proc.time()[[3]] - began)
}

study_signature <- function(root, cfg, versions) {
  hashes <- tools::md5sum(c(file.path(root, "R", "study.R"),
                             file.path(root, "config", "study_config.R")))
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(list(config = cfg, versions = versions, code = unname(hashes)), tmp, version = 2)
  unname(tools::md5sum(tmp))
}

run_study <- function(root = getwd(), worker_id = 1L, n_workers = 1L,
                       cfg = load_config(root)) {
  validate_config(cfg)
  stopifnot(n_workers >= 1, worker_id >= 1, worker_id <= n_workers)
  versions <- check_dependencies()
  signature <- study_signature(root, cfg, versions)
  path <- file.path(root, "results", cfg$run_id)
  for (sub in c("replicates", "logs", "manifests", "locks"))
    dir.create(file.path(path, sub), recursive = TRUE, showWarnings = FALSE)
  existing_manifests <- list.files(file.path(path, "manifests"), "[.]rds$", full.names = TRUE)
  for (f in existing_manifests) {
    if (!identical(readRDS(f)$signature, signature))
      stop("This run_id contains a different configuration/code/package version. Set a NEW run_id.")
  }
  manifest <- list(signature = signature, config = cfg, versions = versions,
                   session = utils::sessionInfo(), worker = worker_id, workers = n_workers,
                   timestamp = Sys.time())
  saveRDS(manifest, file.path(path, "manifests", sprintf("worker_%02d.rds", worker_id)))
  log <- file.path(path, "logs", sprintf("worker_%02d.log", worker_id))
  announce <- function(...) {
    line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
    cat(line, "\n"); cat(line, "\n", file = log, append = TRUE)
    flush.console()
  }
  tasks <- make_tasks(cfg)
  tasks <- tasks[(tasks$replicate - 1L) %% n_workers == worker_id - 1L, , drop = FALSE]
  announce("Run", cfg$run_id, ": worker", worker_id, "of", n_workers, ";", nrow(tasks), "datasets")
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, , drop = FALSE]
    dest <- file.path(path, "replicates", paste0(task$task_id, ".rds"))
    if (file.exists(dest)) {
      old <- readRDS(dest)
      if (!identical(old$signature, signature)) stop("Checkpoint mismatch: ", dest)
      announce("Resume: already complete", task$task_id)
      next
    }
    lock <- file.path(path, "locks", task$task_id)
    if (!dir.create(lock, showWarnings = FALSE)) {
      announce("Locked by another run; skipped", task$task_id)
      next
    }
    tryCatch({
      announce(sprintf("[%d/%d]", i, nrow(tasks)), "Starting", task$task_id)
      result <- run_dataset(cfg, task, signature)
      temp <- paste0(dest, ".tmp-", Sys.getpid())
      saveRDS(result, temp)
      if (!file.rename(temp, dest)) stop("Could not finalize checkpoint: ", dest)
      announce("Finished", task$task_id, sprintf("in %.1f seconds; errors=%d; nonconverged=%d",
        result$elapsed_seconds, sum(!result$status$success),
        sum(result$status$success & !result$status$converged, na.rm = TRUE)))
    }, finally = unlink(lock, recursive = TRUE))
  }
  announce("Worker complete. Run scripts/04_summarize.R after both workers finish.")
  invisible(path)
}
