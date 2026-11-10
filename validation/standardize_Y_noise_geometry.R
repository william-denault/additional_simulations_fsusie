# Inspect the known GP noise covariance after the installed linear transforms.
# Scaling factors come from the observed phenotype, as in susiF.
source("R/study.R")
check_dependencies()
z <- readRDS("validation/standardize_Y_diagnostic.rds")
m <- ncol(z$Y)
basis <- fsusieR::remap_data(diag(m), seq_len(m), verbose = FALSE)$Y
dwt <- get("DWT2", envir = asNamespace("fsusieR"))
rows <- list()
for (i in 1:2) {
  fit <- z$fits[[i]]
  a <- sweep(basis, 2, fit$csd_Y_pos, "/")
  w <- dwt(a, filter.number = fit$filter.number, family = fit$family)
  a <- sweep(cbind(w$D, w$C), 2, fit$csd_Yf, "/")
  # Original E has covariance V; its transformed row E %*% a has a' V a.
  v <- crossprod(a, z$V %*% a)
  noise_var <- diag(v)
  corr <- stats::cov2cor(v)
  rows[[i]] <- data.frame(standardize_Y = fit$standardize_Y,
    min_wavelet_noise_variance = min(noise_var), max_wavelet_noise_variance = max(noise_var),
    noise_variance_ratio = max(noise_var) / min(noise_var),
    max_abs_offdiagonal_correlation = max(abs(corr[row(corr) != col(corr)])))
}
out <- do.call(rbind, rows)
utils::write.csv(out, "validation/standardize_Y_noise_geometry.csv", row.names = FALSE)
print(out, row.names = FALSE)
