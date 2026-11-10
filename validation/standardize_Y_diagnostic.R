# Reproduce the seed-666 GP vignette, holding all settings except the two
# explicitly varied options fixed. This does not change any simulation run.
source("R/study.R")
check_dependencies()
set.seed(666)
n <- 100L
m <- 60L
p <- 12L
maf <- 0.05 + 0.45 * runif(p)
snpids <- paste0("SNP-", seq_len(p))
cpgids <- paste0("CpG-", seq_len(m))
X <- (runif(n * p) < maf) + (runif(n * p) < maf)
X <- matrix(X, n, p, byrow = TRUE)
storage.mode(X) <- "double"
X[, 4] <- X[, 3] + 0.03 * rnorm(n)
colnames(X) <- snpids
F <- matrix(0, p, m, dimnames = list(snpids, cpgids))
F[1, 16:30] <- 2.3
F[9, 16:30] <- -2.3
F[3, 46:60] <- 2
causal <- c(1L, 3L, 9L)
distance <- outer(seq_len(m), seq_len(m), "-")
V <- 0.85 * exp(-0.5 * (distance / 12)^2) + 0.15 * diag(m)
Y <- X %*% F + matrix(rnorm(n * m), n, m) %*% chol(V)
colnames(Y) <- cpgids
settings <- expand.grid(standardize_Y = c(FALSE, TRUE), filter_cs = c(FALSE, TRUE))
rows <- vector("list", nrow(settings))
fits <- vector("list", nrow(settings))
for (i in seq_len(nrow(settings))) {
  fit <- fsusieR::susiF(Y, X, L = 10, filter_cs = settings$filter_cs[i], verbose = FALSE,
    standardize_Y = settings$standardize_Y[i])
  post <- posterior_summary(fit, "fSuSiE", p)
  score <- score_sets(post$cs, causal)
  stopifnot(identical(fit$standardize_Y, settings$standardize_Y[i]))
  rows[[i]] <- cbind(settings[i, ], n_cs = score["n_cs"],
    covered_cs = score["covered_cs"], false_cs = score["n_cs"] - score["covered_cs"],
    cs_power = score["cs_power"], converged = post$converged,
    sets = paste(vapply(post$cs, function(s) paste(s, collapse = ","), character(1)), collapse = "; "))
  fits[[i]] <- fit
}
tab <- do.call(rbind, rows)
rownames(tab) <- NULL
utils::write.csv(tab, "validation/standardize_Y_diagnostic.csv", row.names = FALSE)
saveRDS(list(settings = settings, fits = fits, X = X, Y = Y, F = F, V = V),
  "validation/standardize_Y_diagnostic.rds")
cat("Installed default standardize_Y:", formals(fsusieR::susiF)$standardize_Y, "\n")
print(tab, row.names = FALSE)
# Check that the study wrapper forwards the new option and keeps L = 5.
cfg <- load_config()
dat <- simulate_dataset(cfg, make_tasks(cfg)[1, ])
method <- method_grid(cfg, "iid")
method <- method[method$engine == "fSuSiE", ]
fit <- fit_method(cfg, dat, method)
stopifnot(isTRUE(fit$standardize_Y), method$method == "fSuSiE_standardized", cfg$fit_L == 5L)
cat("PASS: study wrapper forwards standardize_Y = TRUE with configured L = 5.\n")
