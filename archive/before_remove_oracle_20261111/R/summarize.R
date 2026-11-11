# Tables and scientific plots from checkpoints. Does not rerun model fits.

group_rows <- function(x, keys, fun) {
  if (!nrow(x)) return(data.frame())
  ids <- do.call(interaction, c(x[keys], list(drop = TRUE, lex.order = TRUE)))
  pieces <- lapply(split(seq_len(nrow(x)), ids), function(i)
    cbind(x[i[1], keys, drop = FALSE], fun(x[i, , drop = FALSE])))
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

cluster_mean_se <- function(x, cluster) {
  a <- tapply(x, cluster, mean)
  if (length(a) < 2) NA_real_ else sd(a) / sqrt(length(a))
}

ratio_summary <- function(num, den, cluster) {
  a <- as.numeric(tapply(num, cluster, sum))
  b <- as.numeric(tapply(den, cluster, sum))
  if (sum(b) == 0) return(c(value = NA_real_, se = NA_real_))
  r <- sum(a) / sum(b)
  se <- if (length(a) < 2) NA_real_ else sd(a - r * b) / (sqrt(length(a)) * mean(b))
  c(value = r, se = se)
}

aggregate_pips <- function(x, keys) {
  group_rows(x, c(keys, "threshold"), function(a) data.frame(
    n_success = nrow(a), power = mean(a$power),
    power_mcse = cluster_mean_se(a$power, a$replicate),
    fdr = mean(a$fdp), fdr_mcse = cluster_mean_se(a$fdp, a$replicate),
    pooled_fdr = if (sum(a$discoveries)) sum(a$fp) / sum(a$discoveries) else NA_real_,
    discovery_probability = mean(a$discoveries > 0), mean_discoveries = mean(a$discoveries)))
}

aggregate_coverage <- function(x, keys) {
  group_rows(x, c(keys, "nominal"), function(a) {
    # Coverage is conditional on returning at least one CS. Zero-CS fits
    # contribute neither to its point estimate nor its Monte Carlo interval.
    has_cs <- a$n_cs > 0
    r <- ratio_summary(a$covered_cs[has_cs], a$n_cs[has_cs], a$replicate[has_cs])
    data.frame(n_success = nrow(a), n_with_cs = sum(has_cs), cs_coverage = unname(r[1]),
      coverage_mcse = unname(r[2]), total_cs = sum(a$n_cs),
      mean_cs = mean(a$n_cs), cs_availability = mean(a$n_cs > 0),
      mean_cs_size = if (sum(a$n_cs)) sum(a$total_size) / sum(a$n_cs) else NA_real_,
      cs_power = mean(a$cs_power), cs_power_mcse = cluster_mean_se(a$cs_power, a$replicate))
  })
}

compact_metrics <- function(status, pips, coverage, cfg, keys) {
  d <- group_rows(status, keys, function(a) data.frame(
    attempted = nrow(a), successful = sum(a$success), errors = sum(!a$success),
    nonconverged = sum(a$success & !a$converged, na.rm = TRUE),
    mean_fit_seconds = mean(a$fit_seconds),
    mean_prior_seconds = mean(a$prior_seconds)))
  pp <- aggregate_pips(pips[abs(pips$threshold - cfg$table_pip_threshold) < 1e-10, ], keys)
  cc <- aggregate_coverage(coverage[abs(coverage$nominal - cfg$coverage) < 1e-10, ], keys)
  if (nrow(pp)) d <- merge(d, pp, by = keys, all.x = TRUE)
  if (nrow(cc)) {
    cc$n_success <- NULL
    d <- merge(d, cc, by = keys, all.x = TRUE)
  }
  d
}

method_labels <- function(a) {
  labels <- c(diagonal = "mvSuSiE: diagonal", canonical = "mvSuSiE: canonical",
               mash = "mvSuSiE: mash mixture", oracle_shape = "mvSuSiE: oracle shape")
  out <- unname(labels[a$prior])
  out[a$method %in% c("fSuSiE_default", "fSuSiE_standardized")] <- "fSuSiE"
  out[a$method == "fSuSiE_shuffled"] <- "fSuSiE: shuffled (IID)"
  factor(out, levels = c(unname(labels), "fSuSiE", "fSuSiE: shuffled (IID)"))
}

write_reviewer_table <- function(metrics, path, cfg) {
  methods <- unique(metrics[c("method", "prior", "residual")])
  prior_names <- c(diagonal = "Diagonal", canonical = "Canonical",
                     mash = "mash canonical mixture", oracle_shape = "Oracle shapes",
                     functional_default = "Default functional prior")
  residual_names <- c(true_fixed = "True covariance, fixed",
    estimated_empirical_init = "Estimated, empirical initialization",
    estimated_true_init = "Estimated, true initialization",
    empirical_fixed = "Empirical covariance, fixed", fsusie_default = "fSuSiE residual estimation")
  table <- data.frame(
    Method = ifelse(methods$method %in% c("fSuSiE_default", "fSuSiE_standardized"), "fSuSiE",
      ifelse(methods$method == "fSuSiE_shuffled", "fSuSiE, shuffled (IID only)", "mvSuSiE")),
    Prior = unname(prior_names[methods$prior]),
    Residual_covariance = unname(residual_names[methods$residual]))
  for (shape in cfg$effect_types) for (noise in cfg$noise_types) {
    a <- metrics[metrics$effect_type == shape & metrics$noise_type == noise, ]
    a <- a[match(methods$method, a$method), ]
    val <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
    text <- paste(val(a$cs_coverage), val(a$power), val(a$fdr), sep = " / ")
    text[is.na(a$method)] <- "--"
    table[[paste(shape, noise, sep = "_")]] <- text
  }
  utils::write.csv(table, file.path(path, "reviewer_table.csv"), row.names = FALSE)
  table
}

plot_results <- function(coverage, cfg, folder) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  palette <- c("mvSuSiE: diagonal" = "#8C564B", "mvSuSiE: canonical" = "#1F77B4",
    "mvSuSiE: mash mixture" = "#9467BD", "mvSuSiE: oracle shape" = "#2CA02C",
    "fSuSiE" = "#D62728")
  shapes <- stats::setNames(rep(16, length(palette)), names(palette))
  policy_labels <- c(true_fixed = "True residual covariance, fixed",
    estimated_empirical_init = "Residual covariance estimated from empirical initialization",
    estimated_true_init = "Residual covariance estimated from true initialization",
    empirical_fixed = "Empirical residual covariance, fixed")
  # Only native reported CSs at the configured primary level (95% by default).
  primary <- coverage[abs(coverage$nominal - cfg$coverage) < 1e-10 &
    (coverage$method %in% c("fSuSiE_default", "fSuSiE_standardized") |
       coverage$engine == "mvSuSiE"), ]
  shown_labels <- c(names(palette)[1:4],
    unique(as.character(method_labels(primary[primary$engine == "fSuSiE", ]))))
  palette <- palette[shown_labels]
  shapes <- shapes[shown_labels]
  nominal_label <- paste0(format(100 * cfg$coverage, trim = TRUE), "%")
  run_label <- paste0(cfg$n_rep, " replicates per condition",
    if (cfg$n_rep < 50L) " - PILOT" else "")
  base_theme <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F0F3F7"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9),
      plot.margin = ggplot2::margin(12, 16, 12, 12))
  facet_labels <- ggplot2::labeller(effect_type = c(step = "Step effect", smooth = "Smooth effect"),
    noise_type = c(iid = "IID noise", gp = "GP noise"))
  # Small horizontal offsets separate overlapping points; ticks still mark L.
  dodge <- ggplot2::position_dodge(width = 0.26)
  for (policy in cfg$residuals) {
    cc <- primary[primary$residual == policy | primary$engine == "fSuSiE", ]
    cc$label <- factor(as.character(method_labels(cc)), levels = names(palette))
    cc$effect_type <- factor(cc$effect_type, levels = cfg$effect_types)
    cc$noise_type <- factor(cc$noise_type, levels = cfg$noise_types)
    for (metric in c("coverage", "cs_power")) {
      d <- cc
      if (metric == "coverage") {
        # Leave missing values at zero-CS conditions so no point
        # implies observed coverage where none can be assessed.
        d$value <- ifelse(d$total_cs > 0, d$cs_coverage, NA_real_)
        d$se <- ifelse(d$total_cs > 0, d$coverage_mcse, NA_real_)
        title <- paste("Observed coverage of", nominal_label, "credible sets")
        ylabel <- "Observed CS coverage (%)"
        definition <- "Coverage: reported CSs containing a causal SNP / all reported CSs. Zero-CS fits excluded."
      } else {
        d$value <- d$cs_power
        d$se <- d$cs_power_mcse
        title <- paste("Power to capture causal SNPs in", nominal_label, "credible sets")
        ylabel <- "CS power (% of causal SNPs captured)"
        definition <- "Power: mean fraction of true causal SNPs in at least one reported CS. Zero-CS fits count as zero."
      }
      d$lower <- 100 * pmax(0, d$value - 1.96 * d$se)
      d$upper <- 100 * pmin(1, d$value + 1.96 * d$se)
      d$value <- 100 * d$value
      g <- ggplot2::ggplot(d, ggplot2::aes(x = true_L, y = value, color = label,
                                          shape = label, group = label))
      if (metric == "coverage")
        g <- g + ggplot2::geom_hline(yintercept = 100 * cfg$coverage,
                                     linetype = "dashed", color = "grey50")
      g <- g +
        ggplot2::geom_point(position = dodge, size = 2.3, na.rm = TRUE) +
        ggplot2::facet_grid(effect_type ~ noise_type, labeller = facet_labels) +
        ggplot2::scale_color_manual(values = palette, drop = FALSE) +
        ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
        ggplot2::scale_x_continuous(breaks = cfg$true_L) +
        ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 25)) +
        ggplot2::coord_cartesian(xlim = range(cfg$true_L) + c(-0.22, 0.22), ylim = c(0, 100)) +
        ggplot2::labs(title = title, subtitle = paste("mvSuSiE:", policy_labels[[policy]]),
          x = "Number of true causal SNPs", y = ylabel, color = NULL, shape = NULL,
          caption = paste0(run_label, "; fSuSiE reference fits repeated across covariance settings.\n",
            definition, "\nPoints only; Monte Carlo uncertainty remains available in the tables.",
            if (metric == "coverage") paste0(" Dashed guide: ", nominal_label, " coverage.") else "")) +
        base_theme + ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                                    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
      for (ext in c("png", "pdf"))
        ggplot2::ggsave(file.path(folder, paste0(metric, "_by_L_", policy, ".", ext)), g,
                         width = 11, height = 8.3, dpi = 160)
    }
  }
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

summarize_study <- function(root = getwd(), run_id = load_config(root)$run_id, make_plots = TRUE) {
  check_dependencies()
  path <- file.path(root, "results", run_id)
  manifests <- list.files(file.path(path, "manifests"), "[.]rds$", full.names = TRUE)
  if (!length(manifests)) stop("No run manifest found: ", path)
  manifest <- readRDS(manifests[1])
  cfg <- manifest$config
  files <- list.files(file.path(path, "replicates"), "[.]rds$", full.names = TRUE)
  if (!length(files)) stop("No completed datasets yet.")
  results <- lapply(files, readRDS)
  if (!all(vapply(results, function(r) identical(r$signature, manifest$signature), logical(1))))
    stop("Mixed checkpoint signatures; do not combine these results.")
  status <- do.call(rbind, lapply(results, `[[`, "status"))
  ps <- cs <- list()
  k <- 0L
  for (r in results) for (i in seq_len(nrow(r$status))) {
    if (!r$status$success[i]) next
    k <- k + 1L
    meta <- r$status[i, c("effect_type", "noise_type", "true_L", "replicate", "task_id",
                           "method", "prior", "residual", "engine", "converged"), drop = FALSE]
    z <- r$outputs[[i]]
    ps[[k]] <- cbind(meta[rep(1, nrow(z$pips)), ], z$pips)
    cs[[k]] <- cbind(meta[rep(1, nrow(z$coverage)), ], z$coverage)
  }
  if (!k) stop("All fits failed. Inspect checkpoint status/error fields.")
  pips <- do.call(rbind, ps)
  coverage <- do.call(rbind, cs)
  keys <- c("effect_type", "noise_type", "true_L", "method", "prior", "residual", "engine")
  pooled_keys <- setdiff(keys, "true_L")
  pp <- aggregate_pips(pips, keys)
  cc <- aggregate_coverage(coverage, keys)
  detail <- compact_metrics(status, pips, coverage, cfg, keys)
  pooled <- compact_metrics(status, pips, coverage, cfg, pooled_keys)
  tables <- file.path(root, "tables", run_id)
  figures <- file.path(root, "figures", run_id)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  export <- list(fit_diagnostics = status, summary_by_L = detail, summary_pooled_L = pooled,
                 power_fdr_curves = pp, coverage_curves = cc,
                 cs_metrics_by_L = cc[abs(cc$nominal - cfg$coverage) < 1e-10, ])
  for (name in names(export))
    utils::write.csv(export[[name]], file.path(tables, paste0(name, ".csv")), row.names = FALSE)
  # Nonconverged fits remain in the primary results and are never hidden.
  # This companion table checks whether excluding them changes the conclusion.
  good_status <- status[status$success & !is.na(status$converged) & status$converged, ]
  if (nrow(good_status)) {
    converged <- compact_metrics(good_status,
      pips[!is.na(pips$converged) & pips$converged, ],
      coverage[!is.na(coverage$converged) & coverage$converged, ], cfg, keys)
    utils::write.csv(converged, file.path(tables, "summary_converged_only.csv"), row.names = FALSE)
  }
  table <- write_reviewer_table(pooled, tables, cfg)
  expected <- nrow(make_tasks(cfg))
  complete <- nrow(unique(status["task_id"])) == expected
  total_seconds <- sum(vapply(results, `[[`, numeric(1), "elapsed_seconds"))
  expected_main_seconds <- mean(vapply(results, `[[`, numeric(1), "elapsed_seconds")) *
    length(cfg$effect_types) * length(cfg$noise_types) * length(cfg$true_L) * load_config(root)$n_rep
  note <- c(paste("Run:", run_id),
    sprintf("Completed datasets: %d / %d. Replicates per condition: %d.", length(results), expected, cfg$n_rep),
    sprintf("Profile points: %d. GP length scale: %.3g points (fraction %.3g); nugget: %.3g; marginal noise variance: %.3g.",
      cfg$positions, cfg$gp_length_fraction * cfg$positions, cfg$gp_length_fraction,
      cfg$gp_nugget, cfg$noise_sd^2),
    sprintf("Fit errors: %d. Nonconverged fits: %d.", sum(!status$success),
             sum(status$success & !status$converged, na.rm = TRUE)),
    sprintf("Sum of dataset runtimes: %.1f min (not parallel wall time).", total_seconds / 60),
    sprintf("Projected full configured study: %.1f min sequential; %.1f min ideal two-worker time (contention may increase this).",
             expected_main_seconds / 60, expected_main_seconds / 120),
    "Table cells: 95% CS coverage / SNP power / SNP FDR at the configured PIP threshold.",
    paste("PIP threshold:", cfg$table_pip_threshold),
    "The compact table pools L=1,2,3 with equal replicate counts when complete. See summary_by_L.csv.",
    "No CS means undefined coverage, not zero coverage. Read CS availability and size alongside coverage.",
    "Plots: number of true causal SNPs versus observed 95% CS coverage or CS power (both in percent).",
    paste0("Plots include fSuSiE (standardize_Y = ", isTRUE(cfg$fsusie_standardize_Y),
      ") and the four mvSuSiE priors only. Shuffled fits are not plotted."),
    "Coverage excludes zero-CS fits. CS power is the mean fraction of causal SNPs in the union of reported CSs; zero-CS fits count as zero power.",
    "Plotted values and contributing fit counts are exported in cs_metrics_by_L.csv; the PIP-based table is unchanged.",
    "FDR = mean(FP / max(1, discoveries)) across datasets; pooled_fdr is also exported.",
    "Nonconverged fits are included and flagged; failed fits are excluded from metrics and counted.",
    "Intervals and Monte Carlo SEs use replicate clusters. Coverage precision also depends on the number of reported CSs.",
    paste("Packages:", paste(names(manifest$versions), manifest$versions, collapse = "; ")))
  writeLines(note, file.path(tables, "RUN_NOTES.txt"))
  heading <- if (cfg$n_rep < 50L || !complete) "PREVIEW / PILOT RESULTS" else "Simulation summary"
  header <- paste0("<th>", html_escape(names(table)), "</th>", collapse = "")
  rows <- apply(table, 1, function(a) paste0("<tr><td>", paste(html_escape(a), collapse = "</td><td>"), "</td></tr>"))
  html <- c('<!doctype html><html><head><meta charset="utf-8"><title>fSuSiE simulation summary</title>',
    '<style>body{font:15px system-ui;max-width:1500px;margin:40px auto;padding:0 24px;color:#182334}h1{font-size:26px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border:1px solid #d7dee8;padding:9px;text-align:left}th{background:#eaf0f7}tr:nth-child(even){background:#f7f9fc}li{margin:7px 0}</style></head><body>',
    paste0("<h1>", heading, "</h1><ul><li>", paste(html_escape(note), collapse = "</li><li>"), "</li></ul>"),
    paste0("<table><thead><tr>", header, "</tr></thead><tbody>"), rows, "</tbody></table></body></html>")
  writeLines(html, file.path(tables, "reviewer_table.html"))
  if (make_plots) plot_results(cc, cfg, figures)
  cat(paste(note, collapse = "\n"), "\nTables:", tables, "\nFigures:", figures, "\n")
  invisible(list(status = status, pips = pp, coverage = cc, table = table, notes = note))
}
