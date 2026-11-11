library(fsusieR)
library(mvsusieR)
library(mashr)
library(susieR)
data("N3finemapping", package="susieR")

get_sim = function(L=1,
                   f_type=c("smooth","step"),
                   pve=.2,
                   n=100,
                   n_pos=100,
                   noise_type=c("iid","GP"),
                   gp_length_scale=.2*n_pos,
                   gp_nugget=.15){

  f_type = match.arg(f_type)
  noise_type = match.arg(noise_type)

  stopifnot(L >= 1,
            n >= 2,
            n_pos >= 40,
            pve > 0,
            pve < 1,
            gp_length_scale > 0,
            gp_nugget > 0,
            gp_nugget <= 1)

  flist = list()

  X = N3finemapping$X
  stopifnot(n <= nrow(X))

  # Sample individuals from the full dataset.
  X = X[sample.int(nrow(X), n), 1:100, drop=FALSE]
  colnames(X) = paste0("SNP-", 1:ncol(X))

  # Remove constant SNPs, including when there is only one.
  constant = which(apply(X, 2, var) < 1e-7)

  if (length(constant) > 0){
    X = X[, -constant, drop=FALSE]
  }

  stopifnot(L <= ncol(X))

  # Generate effects within randomly positioned 40-point blocks.
  for (l in 1:L){

    start = sample.int(n_pos - 40 + 1, 1)
    block = start:(start + 39)

    effect = rep(0, n_pos)

    if (f_type == "smooth"){
      # Cosine bump: returns smoothly to zero at both boundaries.
      effect[block] =
        1.2 * (1 - cos(seq(0, 2*pi, length.out=40))) / 2
    }

    if (f_type == "step"){
      effect[block] = 2
    }

    flist[[l]] = effect
  }

  F = matrix(0, ncol(X), n_pos)
  true_pos = sample.int(nrow(F), L)

  for (l in 1:L){
    F[true_pos[l], ] = flist[[l]]
  }

  rownames(F) = colnames(X)
  colnames(F) = paste0("CpG-", 1:n_pos)

  Y_true = X %*% F

  # Signal variance averaged across positions.
  signal_var = mean(apply(Y_true, 2, var))
  stopifnot(is.finite(signal_var), signal_var > 0)

  # Target PVE = signal_var / (signal_var + noise_var).
  noise_var = signal_var * (1 - pve) / pve

  if (noise_type == "iid"){
    V = noise_var * diag(n_pos)
  }

  if (noise_type == "GP"){

    distance = outer(1:n_pos, 1:n_pos, "-")
    K = exp(-.5 * (distance / gp_length_scale)^2)

    V = noise_var * (
      (1 - gp_nugget) * K +
        gp_nugget * diag(n_pos)
    )
  }

  # Independent individuals; covariance V across positions.
  E = matrix(rnorm(n*n_pos), n, n_pos) %*% chol(V)

  Y = Y_true + E

  colnames(Y) = colnames(F)
  colnames(E) = colnames(F)
  dimnames(V) = list(colnames(F), colnames(F))

  # Realized signal/noise variance ratio in this sample.
  pve_realized =
    signal_var / (signal_var + mean(apply(E, 2, var)))

  return(list(
    X=X,
    Y=Y,
    F=F,
    Y_true=Y_true,
    E=E,
    true_pos=true_pos,
    V=V,
    pve=pve,
    pve_realized=pve_realized
  ))
}
run_methods = function(sim,
                       L_fit=10,
                       coverage=.95,
                       min_abs_corr=.5,
                       pip_threshold=.9,
                       residual=c("estimated",
                                  "true_fixed",
                                  "estimated_true_init")){

  residual = match.arg(residual)

  X = sim$X
  Y = sim$Y
  p = ncol(X)
  r = ncol(Y)
  truth = sim$true_pos

  # Same covariance settings for the three mvSuSiE priors.
  run_mv = function(prior){

    mvsusieR::mvsusie(
      X, Y,
      L=L_fit,
      prior_variance=prior,
      estimate_prior_variance=TRUE,
      estimate_prior_mixture_weights=TRUE,
      residual_variance=
        if (residual == "estimated") NULL else sim$V,
      estimate_residual_variance=
        (residual != "true_fixed"),
      standardize=TRUE,
      coverage=coverage,
      min_abs_corr=min_abs_corr,
      n_thread=1,
      verbose=FALSE
    )
  }

  fits = list()

  # 1. fSuSiE
  fits[["fSuSiE"]] = fsusieR::susiF(
    Y, X,
    L=L_fit,
    standardize_Y=TRUE,
    filter_cs=TRUE,
    cov_lev=coverage,
    min_purity=min_abs_corr,
    verbose=FALSE
  )

  # 2. mvSuSiE: diagonal prior
  prior_diag = mvsusieR::create_mixture_prior(
    mixture_prior=list(
      matrices=list(.2 * diag(apply(Y, 2, var))),
      weights=1
    )
  )

  fits[["mvSuSiE_diagonal"]] = run_mv(prior_diag)

  # 3. mvSuSiE: canonical prior
  prior_canon = mvsusieR::create_mixture_prior(R=r)

  fits[["mvSuSiE_canonical"]] = run_mv(prior_canon)

  # 4. mvSuSiE: mash prior
  Z = susieR::calc_z(X, Y, center=TRUE, scale=TRUE)

  mash_data = mashr::mash_set_data(
    Bhat=Z,
    Shat=matrix(1, p, r)
  )

  mash_fit = mashr::mash(
    mash_data,
    Ulist=mashr::cov_canonical(mash_data),
    outputlevel=0,
    verbose=FALSE
  )

  prior_mash = mvsusieR::create_mixture_prior(
    fitted_g=mash_fit$fitted_g
  )

  fits[["mvSuSiE_mash"]] = run_mv(prior_mash)

  # Collect metrics and PIPs for this simulation.
  stats = list()
  pips = list()
  sets = list()

  for (method in names(fits)){

    fit = fits[[method]]

    if (method == "fSuSiE"){
      cs = fit$cs
    } else {
      cs = fit$sets$cs
    }

    # Remove empty placeholders and duplicate sets.
    cs = cs[lengths(cs) > 0]

    cs = lapply(cs, function(s){
      sort(unique(as.integer(s)))
    })

    cs_id = vapply(cs, function(s){
      paste(s, collapse=",")
    }, character(1))

    cs = cs[!duplicated(cs_id)]

    pip = as.numeric(fit$pip)

    stopifnot(length(pip) == p,
              all(is.finite(pip)),
              all(pip >= -1e-10 & pip <= 1+1e-10))

    # Does each CS contain at least one true causal SNP?
    hit = vapply(cs, function(s){
      any(s %in% truth)
    }, logical(1))

    selected = which(pip >= pip_threshold)

    stats[[method]] = data.frame(
      method=method,
      residual_policy=
        if (method == "fSuSiE") "fSuSiE" else residual,

      n_cs=length(cs),
      covered_cs=sum(hit),

      # Fraction of reported CSs containing a causal SNP.
      coverage=
        if (length(cs) > 0) mean(hit) else NA_real_,

      # Fraction of true SNPs captured by at least one CS.
      cs_power=mean(truth %in% unlist(cs)),

      pip_threshold=pip_threshold,

      # SNP selection using the PIP threshold.
      pip_power=mean(truth %in% selected),
      pip_fdp=sum(!selected %in% truth) /
        max(1, length(selected)),

      converged=
        if (is.null(fit$converged)) NA else isTRUE(fit$converged)
    )

    pips[[method]] = data.frame(
      method=method,
      snp=colnames(X),
      snp_index=seq_len(p),
      causal=seq_len(p) %in% truth,
      PIP=pip
    )

    sets[[method]] = cs
  }

  return(list(
    summary=do.call(rbind, stats),
    PIP=do.call(rbind, pips),
    CS=sets
  ))
}
