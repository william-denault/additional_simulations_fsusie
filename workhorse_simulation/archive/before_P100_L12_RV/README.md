# fSuSiE / mvSuSiE workhorse simulation

This is a separate project based on `get_sim()` and `run_methods()`. It does not reuse the earlier toy-data results. **No method receives the generating residual covariance or residual variance.**

## Run in RStudio

1. Open `workhorse_simulation.Rproj`.
2. Optionally source `00_check_and_pilot.R` as a Background Job. It runs one replicate per condition and writes a separate, clearly labelled pilot.
3. In the **main RStudio console**, run:

   ```r
   source("05_launch_two_jobs.R")
   ```

   This starts two Background Jobs. Alternatively, open `01_run_worker_1.R` and `02_run_worker_2.R` and choose **Source as Background Job** for each. Do not launch script 05 itself as a Background Job: the RStudio jobs API is available in the main session.
4. After both jobs finish, run:

   ```r
   source("04_summarize_and_plot.R")
   ```

   The summarizer can itself run as a Background Job. It also accepts incomplete runs and labels their figures `PARTIAL RUN`.

`03_run_sequential.R` runs the same study in one process. Use either the two-worker route or the sequential route. Each script locates the project from its own file path, so it also works when sourced with an absolute path from another working directory.

Required packages: `susieR`, `fsusieR`, `mvsusieR`, `mashr`; `rstudioapi` for the launcher. Plotting uses base R. Nothing is installed automatically.

## Study settings

Edit `config.R` before starting. The supplied study has 50 replicates for each of 12 conditions: 1, 2 or 3 causal SNPs; step or smooth effects; IID or GP noise. There are 600 datasets and 2,400 fits. Each of the two workers processes 25 replicates per condition (300 datasets), with one numerical thread requested per worker.

- 100 individuals sampled from all available rows of `susieR::N3finemapping$X`, using its first 200 SNP columns. Monomorphic columns in the sample are removed, while original SNP labels are retained.
- 60 ordered positions. Each causal SNP affects a randomly placed block of 40 positions. Step effects equal 2 within the block. Smooth effects use a raised cosine bump, `1.2 * (1 - cos(t)) / 2`, with zero values at the block endpoints. This avoids introducing step discontinuities into the smooth condition.
- Target PVE is 20%, defined using average across-position signal variance: `signal_var / (signal_var + noise_var)`. Actual Gaussian noise draws are not rescaled to force exactly 20%; realized component PVE is saved.
- Independent noise between individuals. Within each individual, IID noise uses `V = noise_var * I`; GP noise uses `V[j,k] = noise_var * (0.85 * exp(-(j-k)^2 / (2*12^2)) + 0.15 * I[j=k])`. The length scale is **12 positions**, or 20% of the profile; the nugget is 15% of marginal noise variance.
- The same replicate and causal count use matching sampled individuals, causal SNPs and block starts across shape/noise conditions. Different replicates have separate fixed seeds. Worker assignment does not change the generated datasets.

The simulation function retains the true covariance for checking generation. The runner deliberately strips it, along with the true effects and residuals, before calling `run_methods()`. True SNP indices are used only for scoring after fitting. No true causal count is passed as the fitting limit.

## Four fits

All methods fit up to **10 effects**, report 95% credible sets and use a 0.5 minimum absolute genotype-correlation purity threshold. Empty and duplicate CSs are removed consistently.

| Plot label | Effect prior / fitting specification | Residual noise |
| --- | --- | --- |
| fSuSiE | `susiF(..., standardize_Y=TRUE, filter_cs=TRUE)` | fSuSiE's own residual-variance estimation |
| mvSuSiE: diagonal | A single diagonal prior initialized from observed phenotype variances | Estimated from observed data |
| mvSuSiE: canonical | `create_mixture_prior(R=ncol(Y))` | Estimated from observed data |
| mvSuSiE: mash | mash learns weights/scales over canonical covariance patterns using marginal Z-scores | Estimated from observed data |

Every mvSuSiE call explicitly uses `residual_variance=NULL` and `estimate_residual_variance=TRUE`. The three mvSuSiE variants use the same residual estimation procedure, differing in effect-prior specification. Prior variance and mixture weights are estimated in mvSuSiE. The mash variant is an illustration trained on this small locus; it does not learn new spatial covariance matrices or use an external training panel. No oracle prior, true covariance or shuffled fSuSiE is included.

The supplied code preserves the package's handling of 60 positions (including fSuSiE's internal wavelet grid conversion). Standardizing Y is explicit in the fSuSiE call; the displayed method name remains **fSuSiE**.

## Metrics

For each simulation and method:

- **Coverage** = number of reported CSs containing at least one true causal SNP / number of reported CSs. If no CS is reported, coverage is `NA`. The main plot averages these per-simulation fractions over simulations with at least one CS, giving each such simulation equal weight. `pooled_coverage` is also exported but is not used for the main coverage plot.
- **CS power** = number of distinct true causal SNPs contained in the union of all reported CSs / number of true causal SNPs. No CS gives zero power. This measures causal-SNP capture; it does not require one separate CS per causal SNP. Interpret it together with CS size.
- **CS availability** = fraction of successfully fitted simulations reporting at least one CS. Removing zero-CS simulations from coverage must not hide this quantity.
- **PIP power** = fraction of causal SNPs with PIP at least the chosen threshold. **PIP FDP** = false SNP discoveries / `max(1, number of discoveries)`. Average FDP across simulations gives the plotted FDR. Default marked threshold is 0.9. PIP power and CS power are different metrics.

Failed fits are excluded from accuracy metrics and counted. Nonconverged fits remain included and flagged; a converged-only sensitivity table is also written. A failed fit is never treated as a successful fit with zero CSs. All SNP PIPs and original SNP identifiers are retained. R-level warnings are saved; also retain the Background Job output because native numerical messages can bypass R's warning capture. A reported convergence flag does not erase earlier ELBO-decrease notices.

Tables contain Monte Carlo SEs and clipped mean +/- 1.96 SE intervals across replicates. These are approximate intervals, and can be degenerate when all observed values are zero or one. The main plots hide intervals by default to match the SuSiE-mix style; set `show_intervals=TRUE` to display them.

## Outputs

- `results/<run_id>/replicates/`: one checkpoint per dataset, including all four methods' metrics, CSs and PIPs, simulation settings and truth metadata.
- `results/<run_id>/logs/`: worker progress and package/session information.
- `tables/<run_id>/summary_by_condition.csv`: method comparisons by causal count, effect shape, noise type and PVE.
- `tables/<run_id>/metrics_per_simulation.csv`, `PIPs_per_simulation.csv`: detailed measurements.
- `tables/<run_id>/PIP_power_FDR.csv`, `PIP_calibration.csv`: plotted PIP summaries.
- `tables/<run_id>/run_report.txt`: completion, errors, convergence and measured runtime.
- `figures/<run_id>/`: both vector PDF and PNG versions of coverage, CS power, CS availability and CS size versus causal count; separate PIP power-FDR curves for each causal count; PIP calibration for each method, pooling causal counts.

Figures follow the local SuSiE-mix plotting scripts: offset coloured points without connecting lines for coverage/power, light grey grids, PVE headings, scenario labels on the right and a dashed 95% coverage reference. The PIP curves preserve threshold order and mark PIP 0.9. The FDR display defaults to 0-25%; change `fdr_plot_max` to see a wider range. Calibration uses ten equal-width PIP bins and replicate-cluster intervals, accounting for paired causal-count scenarios when pooling them. Source style: `susie_mix/script/sim/plot_simulations.R` and `simulation_plot_helpers.R` in the neighbouring analysis directory.

The pilot has only one replicate per condition and checks execution and output, not method performance. Use its runtime as a rough guide; parallel contention and difficult fits may increase full-run time. With 200 candidate SNPs and 60 positions, the full run can take hours even with two workers; the earlier 12-SNP toy study is not a reliable runtime guide.

## Resume and changes

Rerunning a worker skips its saved datasets. Code hashes, package versions and fitting settings are checked before reuse. **Choose a new `run_id` if you change simulation/fitting settings or packages.** Changing plot display settings or the PIP curve threshold grid does not require refitting.

To restart after closing RStudio, stopping a job or rebooting, open the project and source `05_launch_two_jobs.R` again. Finished datasets are skipped; an interrupted dataset is recomputed. Checkpoints are written to temporary files and renamed only after the complete result is saved, so an interrupted write is not mistaken for a completed dataset. There are no persistent locks to clear. Do not delete saved replicate checkpoints to resume. Stop old workers before launching replacements to avoid duplicate computation.

The small deterministic checks can be run from the project directory with `source("validation/checks.R")`. The background-job pilot additionally checks that all four installed packages fit generated data without access to the true residual covariance.
