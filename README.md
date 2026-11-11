 
## Regenerate all tables and plots

Open `minor_update_fsusie.Rproj`, then run in RStudio:

```r
source("run_get_summary.R")
```
 
To regenerate just one study from the repository root:

```r
# Covariance benchmark only:
source("scripts/04_summarize.R")

# Workhorse only:
source("workhorse_simulation/04_summarize_and_plot.R")
```
 
