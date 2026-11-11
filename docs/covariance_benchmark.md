# Controlled covariance benchmark

The completed run is `main_100_T60_gpL12_stdY`: 100 replicates in each of 12
conditions, for 1,200 datasets. Step and smooth effects are crossed with IID and
GP noise and 1, 2, or 3 causal SNPs. Current reports include 10 fits per dataset:
fSuSiE and three mvSuSiE priors crossed with three residual covariance policies.

## Data generation

- 100 individuals, 12 candidate SNPs, and 60 ordered measurements.
- Binomial genotypes with MAFs between 0.10 and 0.45, centered and scaled.
  SNP-4 is a highly correlated, noncausal proxy for SNP-3.
- Causal SNPs are the first L entries of `(1, 3, 9)`. SNPs 1 and 9 affect the
  second quarter of the profile in opposite directions; SNP 3 affects the fourth quarter.
- Step effects and smooth positive sine bumps have the same support and squared
  norm. Mean squared effect per position is 0.10 per causal SNP. Total genetic
  contribution increases with L; this design does not fix PVE at 20%.
- Marginal noise variance is 1. IID noise has covariance I. GP noise has
  `V[j,k] = 0.85 * exp(-(j-k)^2/(2*12^2)) + 0.15 * I[j=k]`.
  Individuals have independent noise profiles. Length scale is 12 positions.
- Data and fitting limits are matched across methods. All fits use L=5,
  95% credible sets, and minimum absolute genotype-correlation purity 0.5.
  The true causal count is not supplied as the fitting limit.

See `config/study_config.R` and `R/study.R` for the exact implementation.

## Prior and residual covariance factors

| Effect prior | Construction |
|---|---|
| Diagonal | One diagonal covariance based on observed phenotype variances |
| Canonical | `mvsusieR::create_mixture_prior(R=60)` |
| mash | Weights/scales learned from the locus over canonical covariance patterns |

The mash variant does not learn new spatial covariance matrices or use a
genome-wide training panel. Each supplied prior is constructed once per dataset
and reused across residual policies. Mixture weights remain fixed; effect-specific
prior scales are estimated. Genotypes are already standardized, and mvSuSiE uses
`standardize=FALSE` to preserve the supplied prior covariance scale.

| Policy | Initial covariance | Updated during fitting? |
|---|---|---|
| `true_fixed` | True generating residual covariance | No |
| `estimated_empirical_init` | Empirical phenotype covariance, `cov(Y)` | Yes |
| `estimated_true_init` | True generating residual covariance | Yes |

The true-covariance settings are simulation diagnostics. They use information
unavailable in a real analysis. In `fit_method()`, `residual_variance=dat$V`
supplies the true matrix and `estimate_residual_variance` controls updating.
The separate workhorse study uses `residual_variance=NULL` for all mvSuSiE fits.

fSuSiE uses its default functional prior, residual estimation, and filtering,
with `standardize_Y=TRUE`. This standardization is not full GP whitening.
fSuSiE receives the same 60 positions and handles its wavelet-grid conversion
internally. The same fSuSiE fits appear in each covariance-policy figure.

Historical oracle-shape fits remain inside the original checkpoints for
provenance but are excluded from every current derived table and figure.
No shuffled fSuSiE fits are included in the current run.

## Metrics

- **Coverage:** total reported CSs containing at least one true causal SNP divided
  by total reported CSs, pooled across simulations. No sets gives `NA`.
  This differs from the workhorse's mean of per-simulation coverage fractions.
- **CS power:** fraction of distinct causal SNPs in the union of reported sets,
  averaged across simulations. A causal SNP is counted once even if it appears
  in multiple sets. No sets gives zero power.
- **CS availability:** fraction of successful fits that report at least one set.
  Set size and availability should be considered alongside coverage.
- **PIP power:** fraction of true causal SNPs with PIP at or above a threshold.
- **PIP FDR:** mean `FP/max(1, discoveries)` across simulations. `pooled_fdr`
  is exported separately. A noncausal SNP in high LD is still a false discovery.

Failed fits are counted and excluded from accuracy metrics. Nonconverged fits
are included and flagged; `summary_converged_only.csv` provides a sensitivity
analysis. Monte Carlo standard errors account for replicate clusters.

`reviewer_table.csv` pools L=1,2,3 and reports coverage / PIP power / PIP FDR
at PIP 0.9. The coverage and CS-power figures instead separate L on the x-axis.
`cs_metrics_by_L.csv` contains the plotted values and contributing counts.
`coverage_curves.csv` varies nominal coverage for the components selected at
the primary 95% setting, without selecting/filtering components again.

## Reproduction and new fits

From the repository root, run `source("scripts/04_summarize.R")` to regenerate
tables and plots from saved checkpoints. The summarizer checks saved manifests
and excludes historical oracle-shape fits before aggregation.

The simulation entry scripts remain in `scripts/`. For new benchmark fits,
choose a new run ID first: the current prior list differs from the historical
checkpoint configuration. Do not modify saved manifests to bypass this check.
Use the main RStudio console for `scripts/05_launch_two_jobs.R`, or source
workers 01 and 02 as separate Background Jobs. The workhorse has its own launcher.

`tests/check_design_and_metrics.R` checks generation and metrics without fitting
models. `tests/audit_saved_coverage.R` independently checks saved CS membership
against truth, including historical fits retained in the raw checkpoints.
