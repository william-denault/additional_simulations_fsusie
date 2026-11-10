utils::data("N3finemapping", package="susieR")

get_sim = function(L=1,
                   f_type=c("smooth","step"),
                   pve=.2, n=100, n_pos=100,
                   noise_type=c("iid","GP"),
                   gp_length_scale=.2*n_pos, gp_nugget=.15, n_snps=100){

  f_type = match.arg(f_type)
  noise_type = match.arg(noise_type)
  stopifnot(length(L)==1, L>=1, L==as.integer(L),
    length(n)==1, n>=2, n==as.integer(n),
    length(n_pos)==1, n_pos>=40, n_pos==as.integer(n_pos),
    length(pve)==1, pve>0, pve<1,
    gp_length_scale>0, gp_nugget>0, gp_nugget<=1,
    length(n_snps)==1, n_snps>=L, n_snps==as.integer(n_snps))

  flist = list()
  X = N3finemapping$X
  stopifnot(n <= nrow(X), n_snps<=ncol(X))
  X = X[sample.int(nrow(X), n), seq_len(n_snps), drop=FALSE]
  colnames(X) = paste0("SNP-", 1:ncol(X))
  constant = which(apply(X, 2, var) < 1e-7)
  if (length(constant)>0) X = X[, -constant, drop=FALSE]
  stopifnot(L <= ncol(X))

  starts = integer(L)
  for (l in 1:L){
    start = sample.int(n_pos-40+1, 1)
    starts[l] = start
    block = start:(start+39)
    effect = rep(0, n_pos)
    if (f_type=="smooth"){
      effect[block] = 1.2*(1-cos(seq(0, 2*pi, length.out=40)))/2
    }
    if (f_type=="step") effect[block] = 2
    flist[[l]] = effect
  }

  F = matrix(0, ncol(X), n_pos)
  true_pos = sample.int(nrow(F), L)
  for (l in 1:L) F[true_pos[l], ] = flist[[l]]
  rownames(F) = colnames(X)
  colnames(F) = paste0("CpG-", 1:n_pos)
  Y_true = X %*% F

  signal_var = mean(apply(Y_true, 2, var))
  stopifnot(is.finite(signal_var), signal_var>0)
  noise_var = signal_var*(1-pve)/pve
  if (noise_type=="iid") V = noise_var*diag(n_pos)
  if (noise_type=="GP"){
    distance = outer(1:n_pos, 1:n_pos, "-")
    K = exp(-.5*(distance/gp_length_scale)^2)
    V = noise_var*((1-gp_nugget)*K+gp_nugget*diag(n_pos))
  }
  E = matrix(rnorm(n*n_pos), n, n_pos) %*% chol(V)
  Y = Y_true+E
  colnames(Y) = colnames(E) = colnames(F)
  dimnames(V) = list(colnames(F), colnames(F))
  pve_realized = signal_var/(signal_var+mean(apply(E, 2, var)))

  return(list(X=X, Y=Y, F=F, Y_true=Y_true, E=E, true_pos=true_pos,
    V=V, pve=pve, pve_realized=pve_realized, starts=starts,
    signal_var=signal_var, noise_var=noise_var))
}
