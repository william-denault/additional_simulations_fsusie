# The generating covariance is NOT an input to any fit.
capture_fit = function(fun){
  warnings = character()
  start = proc.time()[[3]]
  output = capture.output({
    fit = tryCatch(withCallingHandlers(fun(), warning=function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error=function(e) e)
  })
  # Retain warnings printed through R's output stream too. Native code can
  # also write directly to the job console, bypassing R's warning/output hooks.
  warnings = c(warnings,output[grepl("WARNING|ELBO decreased",output,ignore.case=TRUE)])
  list(fit=fit, success=!inherits(fit,"error"),
       seconds=proc.time()[[3]]-start, warnings=paste(unique(warnings),collapse=" | "))
}

score_fit = function(fit, method, X, truth, pip_threshold, residual_initial=NULL){
  cs = if (method=="fSuSiE") fit$cs else fit$sets$cs
  cs = cs[lengths(cs)>0]
  cs = lapply(cs, function(s) sort(unique(as.integer(s))))
  cs = cs[!duplicated(vapply(cs, function(s) paste(s,collapse=","), character(1)))]
  pip = as.numeric(fit$pip)
  stopifnot(length(pip)==ncol(X), all(is.finite(pip)),
    all(pip>=-1e-10 & pip<=1+1e-10),
    all(unlist(cs)>=1 & unlist(cs)<=ncol(X)))
  hit = vapply(cs, function(s) any(s %in% truth), logical(1))
  selected = which(pip>=pip_threshold)
  residual_change = NA_real_
  if (method!="fSuSiE" && !is.null(residual_initial) && !is.null(fit$sigma2)){
    residual_change = norm(as.matrix(fit$sigma2)-residual_initial,"F")/norm(residual_initial,"F")
  }
  stats = data.frame(n_cs=length(cs), covered_cs=sum(hit),
    total_cs_size=sum(lengths(cs)),
    coverage=if (length(cs)>0) mean(hit) else NA_real_,
    cs_power=mean(truth %in% unlist(cs)),
    mean_cs_size=if (length(cs)>0) mean(lengths(cs)) else NA_real_,
    cs_available=as.numeric(length(cs)>0),
    pip_threshold=pip_threshold, pip_power=mean(truth %in% selected),
    pip_fdp=sum(!selected %in% truth)/max(1,length(selected)),
    residual_change_from_initial=residual_change,
    converged=if (is.null(fit$converged)) NA else isTRUE(fit$converged))
  pips = data.frame(method=method, snp=colnames(X), snp_index=seq_len(ncol(X)),
                    causal=seq_len(ncol(X)) %in% truth, PIP=pip)
  list(stats=stats, PIP=pips, CS=cs)
}

run_methods = function(sim, L_fit=10, coverage=.95,
                       min_abs_corr=.5, pip_threshold=.9){
  X = sim$X
  Y = sim$Y
  p = ncol(X)
  r = ncol(Y)
  truth = sim$true_pos

  # Both variants start at the same data-derived covariance, cov(Y).
  # FALSE holds that initialization fixed; TRUE updates it during fitting.
  # Neither option supplies the true simulation covariance.
  residual_initial = cov(Y)
  run_mv = function(prior,estimate){
    mvsusieR::mvsusie(X, Y, L=L_fit,
      prior_variance=prior, estimate_prior_variance=TRUE,
      estimate_prior_mixture_weights=TRUE,
      residual_variance=NULL, estimate_residual_variance=estimate,
      standardize=TRUE, coverage=coverage, min_abs_corr=min_abs_corr,
      n_thread=1, verbose=FALSE)
  }

  # Build each initial prior once, then give independent copies to the two
  # covariance variants. In particular, mash training is shared within a pair.
  prior_cache = new.env(parent=emptyenv())
  get_prior = function(type){
    if (!exists(type,envir=prior_cache,inherits=FALSE)){
      if (type=="diagonal") prior = mvsusieR::create_mixture_prior(
        mixture_prior=list(matrices=list(.2*diag(apply(Y,2,var))), weights=1))
      if (type=="canonical") prior = mvsusieR::create_mixture_prior(R=r)
      if (type=="mash"){
        Z = susieR::calc_z(X,Y,center=TRUE,scale=TRUE)
        d = mashr::mash_set_data(Bhat=Z,Shat=matrix(1,p,r))
        g = mashr::mash(d,Ulist=mashr::cov_canonical(d),outputlevel=0,verbose=FALSE)
        prior = mvsusieR::create_mixture_prior(fitted_g=g$fitted_g)
      }
      assign(type,prior,envir=prior_cache)
    }
    unserialize(serialize(get(type,envir=prior_cache),NULL))
  }
  methods = list(fSuSiE=function(){
    fsusieR::susiF(Y,X,L=L_fit,standardize_Y=TRUE,filter_cs=TRUE,
      cov_lev=coverage,min_purity=min_abs_corr,verbose=FALSE)
  })
  metadata = list(fSuSiE=data.frame(prior_type="wavelet",residual_policy="fsusie_default",
    estimate_residual_variance=NA))
  for (type in c("diagonal","canonical","mash")) for (estimate in c(FALSE,TRUE)){
    policy = if (estimate) "estimated" else "fixed"
    method = paste("mvSuSiE",type,policy,sep="_")
    methods[[method]] = local({
      prior_type = type; update = estimate
      function() run_mv(get_prior(prior_type),update)
    })
    metadata[[method]] = data.frame(prior_type=type,residual_policy=policy,
      estimate_residual_variance=estimate)
  }

  stats = pips = sets = list()
  for (method in names(methods)){
    result = capture_fit(function(){
      score_fit(methods[[method]](),method,X,truth,pip_threshold,residual_initial)
    })
    if (result$success){
      value = result$fit
      row = value$stats
      pips[[method]] = value$PIP
      sets[[method]] = value$CS
      error = ""
    } else {
      row = data.frame(n_cs=NA_integer_, covered_cs=NA_integer_, total_cs_size=NA_real_,
        coverage=NA_real_, cs_power=NA_real_, mean_cs_size=NA_real_, cs_available=NA_real_,
        pip_threshold=pip_threshold, pip_power=NA_real_, pip_fdp=NA_real_,
        residual_change_from_initial=NA_real_,converged=NA)
      pips[[method]] = data.frame(method=method, snp=colnames(X), snp_index=seq_len(p),
        causal=seq_len(p) %in% truth, PIP=NA_real_)
      sets[method] = list(NULL)
      error = conditionMessage(result$fit)
    }
    stats[[method]] = cbind(data.frame(method=method, success=result$success),metadata[[method]], row,
      seconds=result$seconds, warnings=result$warnings, error=error)
  }
  return(list(summary=do.call(rbind,stats), PIP=do.call(rbind,pips), CS=sets))
}
