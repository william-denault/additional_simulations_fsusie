# Independent audit of saved native 95% CSs. Does not fit models or alter results.
# Rscript tests/audit_saved_coverage.R [project_root] [audit_output_directory]
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], winslash = "/") else getwd()
out_dir <- if (length(args) >= 2L) args[2] else file.path(root, "validation")
config_env <- new.env(parent = baseenv())
sys.source(file.path(root, "config", "study_config.R"), envir = config_env)
run_id <- config_env$study_config$run_id
result_dir <- file.path(root, "results", run_id)
manifest_files <- list.files(file.path(result_dir, "manifests"), "[.]rds$", full.names = TRUE)
stopifnot(length(manifest_files) > 0)
manifest <- readRDS(manifest_files[1])
cfg <- manifest$config
files <- list.files(file.path(result_dir, "replicates"), "[.]rds$", full.names = TRUE)
stopifnot(length(files) > 0)
audit <- vector("list", length(files))
native_sets_checked <- 0L
empty_fits <- 0L

for (k in seq_along(files)) {
  r <- readRDS(files[k])
  stopifnot(identical(r$signature, manifest$signature),
    identical(as.integer(r$causal), c(1L, 3L, 9L)[seq_len(r$task$true_L)]))
  rows <- vector("list", nrow(r$status))
  for (i in seq_len(nrow(r$status))) {
    if (!r$status$success[i]) next
    z <- r$outputs[[i]]
    sets <- z$post$cs
    truth <- r$causal
    nsets <- length(sets)
    stopifnot(all(lengths(sets) > 0L), nrow(z$post$selected_alpha) == nsets)
    if (!nsets) empty_fits <- empty_fits + 1L
    keys <- vapply(sets, function(s) paste(sort(as.integer(s)), collapse = ","), character(1))
    stopifnot(!anyDuplicated(keys))
    hits <- logical(nsets)
    for (j in seq_along(sets)) {
      members <- as.integer(sets[[j]])
      stopifnot(all(members >= 1L & members <= cfg$p), !anyDuplicated(members))
      hits[j] <- length(intersect(members, truth)) > 0L
      a <- as.numeric(z$post$selected_alpha[j, ])
      stopifnot(length(a) == cfg$p, all(is.finite(a)), all(a >= 0),
        abs(sum(a) - 1) < 1e-8, sum(a[members]) >= cfg$coverage - 1e-10)
      # Independently reconstruct the shortest descending-probability prefix.
      ordered <- order(a, decreasing = TRUE)
      cut <- which(cumsum(a[ordered]) >= cfg$coverage - 1e-12)[1]
      stopifnot(!is.na(cut), setequal(members, ordered[seq_len(cut)]))
      native_sets_checked <- native_sets_checked + 1L
    }
    caught <- length(intersect(truth, unique(as.integer(unlist(sets)))))
    saved <- z$coverage[abs(z$coverage$nominal - cfg$coverage) < 1e-10, ]
    stopifnot(nrow(saved) == 1L, saved$n_cs == nsets,
      saved$covered_cs == sum(hits), saved$total_size == sum(lengths(sets)),
      saved$recovered_causal == caught, abs(saved$cs_power - caught / length(truth)) < 1e-12)
    rows[[i]] <- data.frame(effect_type = r$task$effect_type, noise_type = r$task$noise_type,
      true_L = r$task$true_L, replicate = r$task$replicate, method = r$status$method[i],
      n_cs = nsets, covered_cs = sum(hits), recovered_causal = caught,
      coverage_per_fit = if (nsets) mean(hits) else NA_real_, cs_power = caught / length(truth))
  }
  audit[[k]] <- do.call(rbind, rows)
}
raw <- do.call(rbind, audit)
group_keys <- c("effect_type", "noise_type", "true_L", "method")
group_id <- do.call(paste, c(raw[group_keys], sep = "|"))
groups <- split(seq_len(nrow(raw)), group_id)
totals <- do.call(rbind, lapply(groups, function(idx) {
  a <- raw[idx, ]
  nsets <- sum(a$n_cs)
  cbind(a[1, group_keys, drop = FALSE], data.frame(n_fits = nrow(a),
    n_with_cs = sum(a$n_cs > 0), total_cs = nsets, covered_cs = sum(a$covered_cs),
    cs_coverage = if (nsets) sum(a$covered_cs) / nsets else NA_real_,
    mean_fit_coverage = if (nsets) mean(a$coverage_per_fit, na.rm = TRUE) else NA_real_,
    cs_power = mean(a$cs_power)))
}))
rownames(totals) <- NULL
published <- read.csv(file.path(root, "tables", run_id, "cs_metrics_by_L.csv"))
id <- function(x) do.call(paste, c(x[group_keys], sep = "|"))
stopifnot(setequal(id(totals), id(published)))
published <- published[match(id(totals), id(published)), ]
for (field in c("n_with_cs", "total_cs", "cs_coverage", "cs_power"))
  stopifnot(isTRUE(all.equal(totals[[field]], published[[field]], tolerance = 1e-12)))

note <- c(paste("Coverage audit for", run_id),
  sprintf("PASS: %d saved datasets, %d successful fits, %d nonempty native credible sets checked.",
    length(files), nrow(raw), native_sets_checked),
  sprintf("PASS: %d zero-CS fits contribute no coverage denominator and zero CS power.", empty_fits),
  "PASS: saved causal SNP identities agree with the simulation design.",
  "PASS: CS indices are valid, identical sets are deduplicated, and no empty placeholder is counted.",
  "PASS: every saved CS is the shortest descending-alpha prefix reaching the configured 95% level.",
  "PASS: independently recomputed per-fit CS counts, hits, sizes and causal-SNP recovery match saved metrics.",
  sprintf("PASS: all %d method/condition coverage and CS-power values match the plotted summary CSV.", nrow(totals)),
  "Coverage = total reported CSs containing >=1 true causal SNP / total reported CSs (equal weight per CS).",
  "CS power = mean across fits of the fraction of distinct true causal SNPs in the union of reported CSs.",
  "mean_fit_coverage is a separate, equal-weight-per-fit alternative, exported for transparency; it is not plotted.",
  "Audit checks saved posterior summaries and scoring. Full fit objects were not saved; no refitting performed.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(note, file.path(out_dir, paste0("coverage_audit_", run_id, ".txt")))
utils::write.csv(totals, file.path(out_dir, paste0("coverage_audit_", run_id, ".csv")), row.names = FALSE)
cat(paste(note, collapse = "\n"), "\n")
print(totals[totals$method == "fSuSiE_default", ], row.names = FALSE)
