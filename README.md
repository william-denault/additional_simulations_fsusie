# Small fSuSiE / mvSuSiE simulation study

Open the existing `minor_update_fsusie.Rproj` in this directory. The default
study has **100 replicates in each of 12 conditions**: step or smooth effects,
IID or within-person GP noise, and 1, 2, or 3 true causal SNPs. Every method
fits the same upper limit, **L = 5**, so it is not supplied the true causal count.
The current run is `main_100_T60_gpL12_stdY`: **60 profile points**, GP length scale
**12 points (20% of the profile)**, and **1,200 datasets / 12,000 included method fits**.
fSuSiE uses `standardize_Y = TRUE` explicitly and is labeled simply **fSuSiE**.
The shared fitting limit remains L = 5 and default CS filtering is retained.
The completed 32-point, 50-replicate results remain under `main_50`; its previous
configuration is saved in `config/archive/main_50.R`.
The previous 200-replicate configuration is archived in
`config/archive/main_200_T60_gpL12.R`; its checkpoints remain separate.

## Run as two RStudio Background Jobs

From the project console:

```r
# Optional quick check: one replicate in each condition, then preview outputs.
source("scripts/00_check_and_pilot.R")

# Launch the actual 100-replicate study as two independent Background Jobs.
source("scripts/05_launch_two_jobs.R")

# After BOTH jobs finish:
source("scripts/04_summarize.R")
```

The two jobs appear in RStudio's Background Jobs pane. One handles odd replicate
numbers, the other even numbers. They do not copy results into the console's
global environment. Both write to this project, using different checkpoint files.
Keep the project working directory as the jobs' working directory.
If jobs from the previous configuration are still running, stop those jobs
before launching the new run. Changing the configuration does not alter settings
already loaded by a running job.

Alternatively, open `scripts/01_run_worker_1.R` and `scripts/02_run_worker_2.R`
and use RStudio's **Source as Background Job** for each. Do not run the two
scripts one after the other in the same console and expect parallel execution.
`scripts/03_run_sequential.R` is an alternative for one process.

Dependencies: **fsusieR, mvsusieR, susieR, mashr, ggplot2**; the automatic job
launcher also uses **rstudioapi**. The scripts check installed packages and
record versions; they do not install or update packages.

## Files and outputs

```text
config/study_config.R              All design and runtime settings
R/study.R                         Simulation, priors, fitting, scoring, checkpoints
R/summarize.R                     Tables, diagnostics, and plots
scripts/                          Pilot, background jobs, and summary entry points
tests/check_design_and_metrics.R  Checks the design and metric calculations
results/main_100_T60_gpL12_stdY/
  replicates/                     One compact RDS checkpoint per dataset
  manifests/                      Configuration, code signature, package versions
  logs/                           Progress logs, one per worker
  locks/                          Protects against duplicate jobs fitting one dataset
tables/main_100_T60_gpL12_stdY/
  reviewer_table.html             Small readable table for review
  reviewer_table.csv              Same table; cells = coverage / power / FDR
  summary_by_L.csv                Full metrics separately for L = 1, 2, 3
  summary_pooled_L.csv             Metrics pooled across the three true-L settings
  power_fdr_curves.csv            All PIP thresholds, FDR/power and Monte Carlo SEs
  coverage_curves.csv             Nominal versus empirical coverage and MCSEs
  cs_metrics_by_L.csv             Plotted 95% CS coverage/power and contributing counts
  fit_diagnostics.csv             Every fit's error/warnings/convergence/runtime
  summary_converged_only.csv      Sensitivity summary excluding nonconverged fits
  RUN_NOTES.txt                   Completion counts, definitions, timings, versions
figures/main_100_T60_gpL12_stdY/
  coverage_by_L_<residual policy>.png/.pdf
  cs_power_by_L_<residual policy>.png/.pdf
```

Pilot outputs use `<run_id>_pilot` and are prominently labeled
as previews. Each figure has four panels (effect shape x noise), with the number
of true causal SNPs (1, 2, 3) on the x-axis. The y-axis is observed 95% CS coverage
or CS power, expressed as a percentage. Separate figures for each mvSuSiE residual
policy keep the number of curves per panel manageable.
The fSuSiE reference curves are repeated across these figures; these are the
same fits, not separate fSuSiE covariance configurations.
Plots show only fSuSiE and the three mvSuSiE priors. The shuffled fSuSiE
diagnostic remains in the saved results and tables but is excluded from plots.
It is disabled for new runs in the current configuration to avoid unnecessary fits.

Rerun `source("scripts/04_summarize.R")` to regenerate the tables and plots from
saved checkpoints without refitting any models. Earlier `coverage_*` and
`power_fdr_*` curve figures, if present, are retained; the new requested plots
have `_by_L_` in their filenames.

## Matched simulation design

Defaults are **100 individuals, 12 candidate SNPs, and 60 ordered positions**.
The individual and SNP counts follow the small step-effect vignette. This is a small
illustrative study, not a genome-scale benchmark. Set `p = 100` or increase the
number of positions for a larger sensitivity analysis, and change `run_id`.
The position count must be divisible by four for the matched quarter-block
effects. It need not be a power of two: all methods receive the same 60-column
phenotype matrix, and fSuSiE performs its default internal wavelet remapping
to 64 points (verified with the installed version).

Genotypes are binomial with MAFs in [0.10, 0.45], centered and scaled to variance
one. SNP-4 is a highly correlated noncausal allele-count proxy for SNP-3. True
causal SNPs are the first L entries of `(1, 3, 9)`. SNPs 1 and 9 affect the same
block in opposite directions; SNP 3 affects a different block.

The step effects follow the blocks in the supplied vignette: the second and
fourth quarters of the profile. Smooth effects use a positive cosine half-wave
(equivalently a sine bump), rescaled onto the SAME supports. Thus this study
adapts the smooth vignette's shape rather than reusing its different 128-position,
approximately 1,000-SNP dataset. Both shapes have the same sum of squared effects
per causal SNP. The default mean squared effect is 0.10 per position, per causal
SNP; increasing L increases the total genetic contribution.

IID residual covariance is sigma^2 I. GP residual covariance is
`sigma^2 * ((1 - nugget) * exp(-distance^2/(2*length_scale^2)) + nugget*I)`.
The default sigma is 1, nugget is 0.15, and length scale is **12 positions**
(`gp_length_fraction = 0.20` times 60). This is twice the previous relative length
scale (0.10). Correlation is about 0.847 at distance 1, 0.516 at distance 12,
and 0.115 at distance 24. The GP component has variance 0.85 and the independent
component variance 0.15.
Each person has an independent Gaussian vector across positions. There is NO
correlation between individuals. Marginal noise variance is identical in both
noise settings. With `chol(V)` upper triangular, noise is `Z %*% chol(V)`.

All methods see the identical dataset in a replicate. Genotypes and underlying
standard-normal noise draws are also paired across conditions by replicate number.
Neither worker assignment nor execution order changes the data or the seeds.

## Comparisons that address the reviewer

Three effect priors are crossed with three residual-covariance policies:

| Prior | Construction |
|---|---|
| diagonal | Single diagonal effect covariance, as in the prior-specification vignette |
| canonical | Standard canonical covariance mixture |
| mash | Weights/scales fitted to this locus over canonical covariance patterns |

The former oracle-shape comparison has been removed from the active configuration,
derived tables and figures. Original checkpoints retain the historical fits for
reproducibility; the summarizer excludes those fits before aggregation. The mash example learns weights/scales
over supplied canonical patterns, not new spatial covariance matrices. Its
small-locus training sample is a limitation; this is not a benchmark of a
genome-wide learned mvSuSiE prior.

| Residual policy | Initial matrix | Updated during fitting? |
|---|---|---|
| true_fixed | True generating V | No |
| estimated_empirical_init | Package default, cov(Y) | Yes |
| estimated_true_init | True generating V | Yes |

`estimated_true_init` supplies the generating noise covariance as the starting
matrix and sets `estimate_residual_variance=TRUE`, so mvSuSiE can change it during
fitting. This is a simulation-only initialization diagnostic. It is distinct
from the removed oracle-shape effect prior. It is implemented in `fit_method()`
in `R/study.R`: the update flag is selected by the policy, `dat$V` is supplied for
`true_fixed` and `estimated_true_init`, and the argument list is passed to
`mvsusieR::mvsusie`. The separate current `workhorse_simulation` project supplies
`residual_variance=NULL` for all mvSuSiE fits instead.

An optional fourth policy, `empirical_fixed`, fixes the package's initial cov(Y).
Add it to `residuals` to include the complete initialization/fixing comparison
from the smooth vignette.

At fixed V, comparing priors evaluates the effect-prior family. At fixed prior,
comparing residual policies evaluates covariance estimation and initialization.
Each prior is built once per dataset and reused across residual policies. For
this controlled study, mvSuSiE's mixture weights are held fixed to their supplied
values; the per-effect overall prior scales remain estimated. This differs from
the package's default mixture-weight updating. The switch is exposed in config
if a separate sensitivity analysis with weight updating is desired.

X is standardized before simulation. `standardize = FALSE` in mvSuSiE prevents
its additional phenotype-SD rescaling from altering the supplied effect covariance matrices;
intercept fitting remains on. All three residual policies use identical settings.
See the [mvSuSiE function documentation](https://stephenslab.github.io/mvsusieR/reference/mvsusie.html)
and [prior examples](https://stephenslab.github.io/mvsusieR/articles/prior_specification.html).

The **fSuSiE baseline keeps the package's default prior, residual estimation,
CS filtering, optimization, and effect postprocessing**, with
**`standardize_Y = TRUE` explicitly enabled**. The shared fitting L, output
verbosity, and explicitly matched 95% CS / 0.5 purity settings are also set.
This option standardizes the phenotype positions and then the wavelet-coefficient
columns. It changes relative weighting across coefficients and is part of the
fitting specification, rather than a display transformation. It does not estimate
or whiten a full GP residual covariance. The plot/table label remains "fSuSiE";
the configuration, run notes and internal method ID record standardization.
It is a practical baseline, not a cell in the mvSuSiE prior-by-covariance factorial.
Between-method differences cannot be attributed to one factor from that comparison
alone.

There is also an optional **fSuSiE position-shuffling diagnostic, under IID noise
only**. A common random column permutation destroys adjacency while preserving
the genotype design, effect norms, and IID residual distribution. It probes use
of ordering while holding the method fixed. It is deliberately omitted under
GP noise, where shuffling would also change the residual covariance structure.

## Metrics and interpretation

- **SNP power:** fraction of true causal SNPs whose PIP exceeds each threshold.
- **SNP FDR:** average `FP / max(1, number selected)` across replicate datasets.
  No selections give FDP = 0, following the usual definition. The ratio of total
  false selections to total selections is separately exported as `pooled_fdr`.
- **95% CS coverage:** fraction of returned credible sets containing at least one
  true causal SNP. This is conditional on a CS being returned. No CS gives NA
  coverage; it is not evidence of good calibration. Read `cs_availability`,
  `mean_cs`, `mean_cs_size`, and `cs_power` alongside coverage.
  The coverage plots omit conditions with no reported CSs across all replicates;
  zero-CS replicates are also excluded from coverage intervals. `n_with_cs` gives
  the number of fits contributing sets. Coverage pools sets across those fits,
  so each reported set contributes equally.
- **Coverage curves (CSV):** select the reported components at the primary 95% setting
  once, then change CS membership using each component's posterior probabilities.
  Components are NOT reselected or purity-filtered again at every nominal level.
  The 95% point uses the native reported CSs. These curves assess calibration
  conditional on that reporting procedure, not unconditional coverage of every
  simulated effect. This avoids refitting at every nominal level and keeps the
  computation small.
- **CS power:** fraction of true causal SNPs appearing in the union of returned
  95% CSs, averaged across replicates and plotted in percent versus the true
  number of causal SNPs. Each true SNP counts once even if it appears in several
  CSs. A fit returning no CSs has zero power and remains in this average.
  This differs from SNP power at a PIP threshold in the reviewer table.
- **Truth is exact:** a noncausal SNP in high LD with a causal SNP is still a
  false SNP discovery. A CS containing both is counted as covered.

The small reviewer table pools L=1,2,3 with equal numbers of replicates when the
study is complete, and uses PIP >= 0.90. The separate-L results and all threshold
curves are also provided. FDR/power Monte Carlo SEs and approximate coverage SEs
use replicate clusters, respecting pairing when pooling over L. The current run
uses 100 replicates per condition; coverage can still be imprecise when few sets
are returned.

Failed fits are excluded from accuracy metrics but counted in the tables.
Nonconverged fits remain in the main results, are explicitly counted, and have a
separate converged-only sensitivity table. Warnings and errors are saved. Do not
interpret an incomplete run or a method with many failed fits as a complete
comparison. Prior-building time is recorded separately and reused across the
three covariance fits; sum dataset runtimes to estimate total compute time.

## Resume and change settings

Rerunning the same background launcher resumes completed dataset checkpoints.
Each fit's seed is deterministic. If a job is stopped halfway through a dataset,
only that dataset needs repeating. Error results within a finished dataset are
recorded rather than silently retried.

Do not edit configuration or update packages while jobs run. A configuration,
code, or package-version mismatch is rejected on resume; use a new `run_id` for
another analysis. The summarizer reads the saved configuration rather than
silently interpreting old results using new settings. It can produce a clearly
labeled partial preview while jobs are running.

After an abrupt process kill, an empty directory can remain in `results/<run>/locks`.
Only after confirming no job is using that run, remove the corresponding stale
lock directory and rerun the launcher. Normal completion/interruption cleans it up.

To check the design and metric calculations without fitting models:

```r
source("tests/check_design_and_metrics.R")
```
