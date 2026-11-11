# fSuSiE and mvSuSiE: reviewer-response simulations

This repository retains **two completed studies**, their simulation and plotting
code, and the saved checkpoints needed to regenerate the results. No new model
fits are needed to update the figures.

| Study | Replicates per condition | Causal SNPs | Candidate SNPs | Residual covariance comparison |
|---|---:|---|---:|---|
| Covariance benchmark (repository root) | 100 | 1, 2, 3 | 12 | True covariance fixed; estimated from empirical initialization; estimated from true initialization |
| Workhorse (`workhorse_simulation/`) | 250 | 1, 2 | 100 before removal of monomorphic columns | Empirical initialization fixed versus updated; no generating covariance supplied |

Both studies use 100 individuals, 60 ordered measurements, step or smooth effects,
IID or within-individual GP noise, and fSuSiE with `standardize_Y=TRUE`.
The workhorse targets 20% PVE. The covariance benchmark instead fixes per-SNP
effect strength; its total PVE is not held at 20% as the causal count changes.

## Regenerate all tables and plots

Open `minor_update_fsusie.Rproj`, then run in RStudio:

```r
source("run_get_summary.R")
```

This reads the **1,200 benchmark and 2,000 workhorse checkpoints**, regenerates
both studies' tables and figures, and restores the original working directory.
It does not run or restart simulations. The script also works when sourced by
absolute path, or as an RStudio Background Job with this project as working directory.

To regenerate just one study from the repository root:

```r
# Covariance benchmark only:
source("scripts/04_summarize.R")

# Workhorse only:
source("workhorse_simulation/04_summarize_and_plot.R")
```

Required packages: `fsusieR`, `mvsusieR`, `susieR`, `mashr`, and `ggplot2`.
Simulation job launching additionally uses `rstudioapi`. Packages are not
installed or upgraded automatically. The saved run reports record the versions used.

## Where to find the results

| Contents | Covariance benchmark | Workhorse |
|---|---|---|
| Figures (PNG and PDF) | `figures/main_100_T60_gpL12_stdY/` | `workhorse_simulation/figures/workhorse_50_P100_L12_T60_RV/` |
| Tables and run report | `tables/main_100_T60_gpL12_stdY/` | `workhorse_simulation/tables/workhorse_50_P100_L12_T60_RV/` |
| Saved checkpoints | `results/main_100_T60_gpL12_stdY/` | `workhorse_simulation/results/workhorse_50_P100_L12_T60_RV/` |

The `50` in the workhorse run ID is historical: the completed study has
**250 replicates per condition**. Keeping the ID preserves checkpoint compatibility.

Start with the workhorse `coverage_by_L`, `cs_power_by_L`, and
`cs_available_by_L` plots. The benchmark provides separate coverage and power
plots for each residual covariance policy. All method points are filled circles
without connecting lines. Oracle-shape and shuffled fSuSiE comparisons are not shown.

The main coverage summaries differ deliberately: the benchmark pools reported
credible sets; the workhorse averages coverage within each simulation and then
across simulations with at least one set. **No sets means undefined coverage**,
so those conditions have no coverage point. CS power counts no-set fits as zero.
Read coverage together with availability, size, and power.

## Code and study details

- `R/`, `config/`, `scripts/`: covariance benchmark implementation and entry scripts.
- `docs/covariance_benchmark.md`: benchmark design, covariance settings, and scoring.
- `tests/`: benchmark design and saved-CS auditing checks.
- `workhorse_simulation/README.md`: workhorse design, fitting settings, metrics, and restart instructions.
- `workhorse_simulation/R/`: generator, seven method fits, checkpointing, summaries, and plots.
- `workhorse_simulation/validation/`: current validation scripts and extension report.

## Run or extend the workhorse yourself

Open `workhorse_simulation/workhorse_simulation.Rproj` and, in the **main RStudio
console**, run `source("05_launch_two_jobs.R")`. Alternatively, from this
repository's project console use `source("run_simulations.R")`. Finished datasets
are skipped. Do not launch these launcher scripts as Background Jobs themselves.

To extend the run, increase only `n_rep` in `workhorse_simulation/config.R` and
keep its run ID. Use a new run ID for scientific setting or fitting-code changes.
The earlier benchmark is complete; use a **new run ID** for new benchmark fits
because its current configuration excludes historical oracle-shape fits.

## Cleanup archive

Pilots, superseded studies, old source backups, and obsolete diagnostics were
moved outside the repository on 11 November 2026 to:

`../additional_simulations_fsusie_archive/cleanup_2026-11-11/`

The archive contains a file manifest and the original versions of edited files.
Both completed studies' checkpoints remain in their original locations.
Git history and local RStudio session files are preserved. RStudio state is
excluded from version control. Current checkpoint files, compact tables, and
figures are eligible for version control; the large, reproducible
`PIPs_per_simulation.csv` export remains local and is regenerated by the summary script.
