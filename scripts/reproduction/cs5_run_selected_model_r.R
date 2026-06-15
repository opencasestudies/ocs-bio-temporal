#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(zellkonverter)
  library(CoGAPS)
  library(SummarizedExperiment)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

IFN_PROGRAM <- c(
  "IFITM1", "IFI6", "ISG15", "MX1", "IFIT1", "IFIT3",
  "IFI44L", "ISG20", "LY6E", "TRIM22", "OAS1", "OASL"
)

pattern_columns <- function(x) {
  pats <- grep("^Pattern_?[0-9]+$", x, value = TRUE)
  pats[order(as.integer(sub("^Pattern_?", "", pats)))]
}

canonicalize_pattern_names <- function(x) {
  sub("^Pattern_?([0-9]+)$", "Pattern\\1", x)
}

write_json_file <- function(path, payload) {
  write_json(payload, path = path, auto_unbox = TRUE, pretty = TRUE, na = "null")
}

ensure_pattern_matrix_orientation <- function(mat, expected_rows, expected_names) {
  if (is.null(dim(mat))) {
    stop("Expected a 2D matrix from CoGAPS.")
  }
  if (nrow(mat) != expected_rows && ncol(mat) == expected_rows) {
    mat <- t(mat)
  }
  if (nrow(mat) != expected_rows) {
    stop(sprintf("Unexpected matrix shape. Expected %d rows, got %d.", expected_rows, nrow(mat)))
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- expected_names
  }
  mat
}

clean_matrix <- function(mat) {
  mat[is.na(mat)] <- 0
  mat[is.infinite(mat)] <- 0
  mat[mat < 0] <- 0
  mat
}

eta_squared <- function(values, groups) {
  values <- as.numeric(values)
  groups <- as.character(groups)
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- groups[keep]
  if (length(values) == 0) {
    return(NA_real_)
  }
  grand_mean <- mean(values)
  ss_total <- sum((values - grand_mean)^2)
  if (!is.finite(ss_total) || ss_total <= 0) {
    return(0)
  }
  ss_between <- 0
  for (group in unique(groups)) {
    vals <- values[groups == group]
    if (length(vals) > 0) {
      ss_between <- ss_between + length(vals) * (mean(vals) - grand_mean)^2
    }
  }
  as.numeric(ss_between / ss_total)
}

classify_pattern <- function(eta_time, eta_cell) {
  if (!is.finite(eta_time) || !is.finite(eta_cell)) {
    return("unclassified")
  }
  if (eta_cell > 1.5 * eta_time) {
    return("identity-like")
  }
  if (eta_time > 1.5 * eta_cell) {
    return("activity-like")
  }
  "mixed"
}

top_genes_table <- function(feature_loadings, top_n) {
  rows <- list()
  row_idx <- 1
  for (pattern in colnames(feature_loadings)) {
    values <- feature_loadings[, pattern]
    ranked <- order(values, decreasing = TRUE)
    ranked <- ranked[seq_len(min(top_n, length(ranked)))]
    for (rank_idx in seq_along(ranked)) {
      rows[[row_idx]] <- data.frame(
        pattern = pattern,
        rank = rank_idx,
        gene = rownames(feature_loadings)[ranked[rank_idx]],
        weight = as.numeric(values[ranked[rank_idx]]),
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }
  do.call(rbind, rows)
}

safe_spearman <- function(x, y) {
  keep <- is.finite(as.numeric(x)) & is.finite(as.numeric(y))
  if (sum(keep) < 3) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  out <- suppressWarnings(cor.test(as.numeric(x)[keep], as.numeric(y)[keep], method = "spearman", exact = FALSE))
  c(rho = unname(out$estimate), p = out$p.value)
}

summarize_patterns <- function(cell_meta_patterns, top_gene_df) {
  pattern_names <- pattern_columns(colnames(cell_meta_patterns))
  time_order <- unique(
    cell_meta_patterns[order(cell_meta_patterns$day_merged_numeric), c("timepoint_merged", "day_merged_numeric")]
  )$timepoint_merged
  rows <- list()
  for (pattern in pattern_names) {
    groups <- lapply(time_order, function(tp) {
      vals <- as.numeric(cell_meta_patterns[cell_meta_patterns$timepoint_merged == tp, pattern])
      vals[is.finite(vals)]
    })
    groups <- groups[vapply(groups, length, integer(1)) > 1]
    kw_stat <- NA_real_
    kw_p <- NA_real_
    if (length(groups) >= 2) {
      kw <- tryCatch(kruskal.test(groups), error = function(e) NULL)
      if (!is.null(kw)) {
        kw_stat <- unname(kw$statistic)
        kw_p <- kw$p.value
      }
    }
    eta_time <- eta_squared(cell_meta_patterns[[pattern]], cell_meta_patterns$timepoint_merged)
    eta_cell <- eta_squared(cell_meta_patterns[[pattern]], cell_meta_patterns$broad_cell_type)
    if ("ifn_program_score" %in% colnames(cell_meta_patterns)) {
      ifn <- safe_spearman(cell_meta_patterns[[pattern]], cell_meta_patterns$ifn_program_score)
      rho_ifn <- ifn[["rho"]]
      p_ifn <- ifn[["p"]]
    } else {
      rho_ifn <- NA_real_
      p_ifn <- NA_real_
    }
    means_by_time <- aggregate(
      cell_meta_patterns[[pattern]],
      by = list(
        timepoint_merged = cell_meta_patterns$timepoint_merged,
        day_merged_numeric = cell_meta_patterns$day_merged_numeric
      ),
      FUN = mean,
      na.rm = TRUE
    )
    means_by_time <- means_by_time[order(means_by_time$day_merged_numeric), ]
    peak_timepoint <- means_by_time$timepoint_merged[which.max(means_by_time$x)]
    top_genes <- head(top_gene_df$gene[top_gene_df$pattern == pattern], 15)
    ifn_overlap <- sum(top_genes %in% IFN_PROGRAM)
    rows[[pattern]] <- data.frame(
      pattern = pattern,
      kruskal_time_stat = kw_stat,
      kruskal_time_p = kw_p,
      eta_timepoint = eta_time,
      eta_broad_cell_type = eta_cell,
      pattern_class = classify_pattern(eta_time, eta_cell),
      spearman_ifn_score_rho = rho_ifn,
      spearman_ifn_score_p = p_ifn,
      peak_timepoint = as.character(peak_timepoint),
      ifn_top_gene_overlap_top15 = ifn_overlap,
      candidate_ifn_pattern = isTRUE(ifn_overlap >= 3 || (is.finite(rho_ifn) && rho_ifn > 0.35)),
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary$kruskal_time_p_adj <- p.adjust(summary$kruskal_time_p, method = "BH")
  summary$ifn_score_p_adj <- p.adjust(summary$spearman_ifn_score_p, method = "BH")
  summary <- summary[order(summary$candidate_ifn_pattern, summary$eta_timepoint, decreasing = TRUE), ]
  rownames(summary) <- NULL
  summary
}

pattern_means <- function(cell_meta_patterns, pattern_names, group_cols) {
  aggregate(cell_meta_patterns[, pattern_names, drop = FALSE], by = cell_meta_patterns[, group_cols, drop = FALSE], FUN = mean, na.rm = TRUE)
}

effective_pattern_summary <- function(feature_loadings, sample_factors) {
  rows <- list()
  n_effective <- 0
  degenerate <- c()
  for (pattern in colnames(feature_loadings)) {
    gene_vec <- as.numeric(feature_loadings[, pattern])
    cell_vec <- as.numeric(sample_factors[, pattern])
    gene_vec[!is.finite(gene_vec)] <- 0
    cell_vec[!is.finite(cell_vec)] <- 0
    gene_max <- max(gene_vec, na.rm = TRUE)
    cell_max <- max(cell_vec, na.rm = TRUE)
    nonzero_gene_count <- sum(gene_vec > 0)
    effective <- isTRUE(gene_max > 0 && cell_max > 0 && sd(cell_vec) > 1e-8 && nonzero_gene_count >= 10)
    if (effective) {
      n_effective <- n_effective + 1
    } else {
      degenerate <- c(degenerate, pattern)
    }
    rows[[pattern]] <- list(
      pattern = pattern,
      effective = effective,
      gene_loading_max = gene_max,
      cell_score_max = cell_max,
      nonzero_gene_count = nonzero_gene_count
    )
  }
  list(rows = unname(rows), n_effective = n_effective, degenerate = degenerate)
}

redundancy_summary <- function(feature_loadings, top_gene_df) {
  corr <- suppressWarnings(cor(feature_loadings, method = "spearman", use = "pairwise.complete.obs"))
  corr[!is.finite(corr)] <- 0
  if (ncol(corr) <= 1) {
    tri <- numeric()
  } else {
    tri <- corr[upper.tri(corr)]
    tri <- pmax(tri, 0)
  }
  jaccards <- c()
  pattern_names <- colnames(feature_loadings)
  for (i in seq_along(pattern_names)) {
    if (i == length(pattern_names)) {
      next
    }
    for (j in (i + 1):length(pattern_names)) {
      left <- top_gene_df$gene[top_gene_df$pattern == pattern_names[[i]]]
      right <- top_gene_df$gene[top_gene_df$pattern == pattern_names[[j]]]
      union <- union(left, right)
      jaccards <- c(jaccards, if (length(union) == 0) 0 else length(intersect(left, right)) / length(union))
    }
  }
  list(
    within_run_pattern_redundancy_mean = if (length(tri)) mean(tri) else 0,
    within_run_pattern_redundancy_max = if (length(tri)) max(tri) else 0,
    within_run_top_gene_jaccard_mean = if (length(jaccards)) mean(jaccards) else 0,
    within_run_top_gene_jaccard_max = if (length(jaccards)) max(jaccards) else 0
  )
}

trace_table <- function(md) {
  chisq <- md$chisq %||% numeric()
  atoms_a <- md$atomsA %||% numeric()
  atoms_p <- md$atomsP %||% numeric()
  n <- max(length(chisq), length(atoms_a), length(atoms_p), 0)
  if (n == 0) {
    return(data.frame(trace_index = integer(), chisq = numeric(), atomsA = numeric(), atomsP = numeric()))
  }
  data.frame(
    trace_index = seq_len(n),
    chisq = c(chisq, rep(NA_real_, n - length(chisq))),
    atomsA = c(atoms_a, rep(NA_real_, n - length(atoms_a))),
    atomsP = c(atoms_p, rep(NA_real_, n - length(atoms_p)))
  )
}

safe_dim <- function(x) {
  dims <- dim(x)
  if (is.null(dims)) {
    return(NULL)
  }
  as.integer(dims)
}

table_count <- function(tab, label) {
  if (label %in% names(tab)) {
    return(unname(as.integer(tab[[label]])))
  }
  0L
}

option_list <- list(
  make_option("--discovery-h5ad", type = "character", dest = "discovery_h5ad", default = "/home/rstudio/project/data/external/reproduction_archive/cache/gse154386_experimental_discovery_cells_x_genes.h5ad", help = "Cells x genes discovery AnnData for R CoGAPS"),
  make_option("--outdir", type = "character", default = "/home/rstudio/project/GSE154386/cogaps_r_revised_model_K10_seed2_iter2000", help = "Output directory"),
  make_option("--k", type = "integer", default = 10, help = "Selected CoGAPS K"),
  make_option("--seed", type = "integer", default = 2, help = "Selected CoGAPS seed"),
  make_option("--n-iter", type = "integer", dest = "n_iter", default = 2000, help = "Selected CoGAPS iteration count"),
  make_option("--top-genes", type = "integer", dest = "top_genes", default = 50, help = "Top genes exported per pattern"),
  make_option("--cogaps-threads", type = "integer", dest = "cogaps_threads", default = 4, help = "Threads passed to CoGAPS"),
  make_option("--blas-threads", type = "integer", dest = "blas_threads", default = 1, help = "Threads exported to BLAS-like libraries"),
  make_option("--output-frequency", type = "integer", dest = "output_frequency", default = 100, help = "Trace output frequency"),
  make_option("--checkpoint-interval", type = "integer", dest = "checkpoint_interval", default = 500, help = "Checkpoint interval"),
  make_option("--checkpoint-outfile", type = "character", dest = "checkpoint_outfile", default = "", help = "Optional checkpoint output path"),
  make_option("--n-snapshots", type = "integer", dest = "n_snapshots", default = 10, help = "Snapshots per enabled phase"),
  make_option("--snapshot-phase", type = "character", dest = "snapshot_phase", default = "all", help = "Snapshot phase: sampling, equilibration, or all"),
  make_option("--take-pump-samples", action = "store_true", dest = "take_pump_samples", default = TRUE, help = "Collect pump diagnostics"),
  make_option("--no-pump-samples", action = "store_false", dest = "take_pump_samples", help = "Disable pump diagnostics"),
  make_option("--asynchronous-updates", action = "store_true", dest = "asynchronous_updates", default = TRUE, help = "Enable asynchronous updates"),
  make_option("--sync-updates", action = "store_false", dest = "asynchronous_updates", help = "Disable asynchronous updates"),
  make_option("--use-sparse-opt", action = "store_true", dest = "use_sparse_opt", default = TRUE, help = "Enable sparseOptimization"),
  make_option("--no-sparse-opt", action = "store_false", dest = "use_sparse_opt", help = "Disable sparseOptimization"),
  make_option("--force-rerun", action = "store_true", default = FALSE, help = "Rerun even if status=ok exists")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (opt$k <= 0 || opt$n_iter <= 0 || opt$top_genes <= 0) {
  stop("--k, --n-iter, and --top-genes must be positive.")
}

outdir <- normalizePath(opt$outdir, winslash = "/", mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
tag <- sprintf("K%d_seed%d_iter%d", opt$k, opt$seed, opt$n_iter)
stem <- sprintf("cogaps_%s", tag)
result_path <- file.path(outdir, sprintf("%s.rds", stem))
metrics_path <- file.path(outdir, sprintf("%s.metrics.json", stem))
diagnostics_path <- file.path(outdir, sprintf("%s.diagnostics.json", stem))
trace_path <- file.path(outdir, sprintf("%s.trace.csv", stem))
log_path <- file.path(outdir, sprintf("%s.log", stem))
checkpoint_path <- opt$checkpoint_outfile
if (identical(checkpoint_path, "")) {
  checkpoint_path <- file.path(outdir, sprintf("%s.checkpoint.out", stem))
}

if (file.exists(result_path) && file.exists(metrics_path) && !isTRUE(opt$force_rerun)) {
  prior <- tryCatch(fromJSON(metrics_path), error = function(e) NULL)
  if (!is.null(prior) && identical(prior$status, "ok")) {
    message(sprintf("[CACHE] status=ok; skipping %s", tag))
    quit(save = "no", status = 0)
  }
}

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

Sys.setenv(
  OMP_NUM_THREADS = as.character(opt$cogaps_threads),
  OPENBLAS_NUM_THREADS = as.character(opt$blas_threads),
  MKL_NUM_THREADS = as.character(opt$blas_threads),
  VECLIB_MAXIMUM_THREADS = as.character(opt$blas_threads),
  NUMEXPR_NUM_THREADS = as.character(opt$blas_threads)
)

message(sprintf("[RUN] Starting selected R CoGAPS model %s", tag))
message(sprintf("[RUN] discovery_h5ad=%s", opt$discovery_h5ad))
message(sprintf("[RUN] outdir=%s", outdir))
t0 <- Sys.time()

metrics <- list(
  created_at_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
  status = "ok",
  language = "R",
  package = "CoGAPS",
  K = opt$k,
  seed = opt$seed,
  n_iter = opt$n_iter,
  distributed_mode = FALSE,
  sparseOptimization = isTRUE(opt$use_sparse_opt),
  nThreads = opt$cogaps_threads,
  blas_threads = opt$blas_threads,
  outputFrequency = opt$output_frequency,
  checkpointInterval = opt$checkpoint_interval,
  checkpointOutFile = checkpoint_path,
  nSnapshots = opt$n_snapshots,
  snapshotPhase = opt$snapshot_phase,
  takePumpSamples = isTRUE(opt$take_pump_samples),
  asynchronousUpdates = isTRUE(opt$asynchronous_updates),
  discovery_h5ad = opt$discovery_h5ad,
  result_path = result_path,
  metrics_path = metrics_path,
  diagnostics_path = diagnostics_path,
  trace_path = trace_path,
  log_path = log_path
)

tryCatch({
  sce <- zellkonverter::readH5AD(opt$discovery_h5ad, reader = "R")
  data_matrix <- SummarizedExperiment::assay(sce, "X")
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "double"
  data_matrix[is.na(data_matrix)] <- 0
  data_matrix[is.infinite(data_matrix)] <- 0
  data_matrix[data_matrix < 0] <- 0

  gene_names <- rownames(sce)
  cell_names <- colnames(sce)
  if (is.null(gene_names)) {
    gene_names <- sprintf("gene_%d", seq_len(nrow(sce)))
    rownames(sce) <- gene_names
  }
  if (is.null(cell_names)) {
    cell_names <- sprintf("cell_%d", seq_len(ncol(sce)))
    colnames(sce) <- cell_names
  }
  rownames(data_matrix) <- gene_names
  colnames(data_matrix) <- cell_names

  openmp_support <- compiledWithOpenMPSupport()
  if (opt$cogaps_threads > 1 && !openmp_support) {
    stop("The installed CoGAPS build does not report OpenMP support; refusing a multithreaded benchmark.")
  }

  params <- CogapsParams(
    nPatterns = opt$k,
    nIterations = opt$n_iter,
    seed = opt$seed,
    sparseOptimization = isTRUE(opt$use_sparse_opt)
  )
  params@takePumpSamples <- isTRUE(opt$take_pump_samples)
  set.seed(opt$seed)

  result <- CoGAPS(
    data_matrix,
    params,
    nThreads = opt$cogaps_threads,
    messages = TRUE,
    outputFrequency = opt$output_frequency,
    checkpointOutFile = checkpoint_path,
    checkpointInterval = opt$checkpoint_interval,
    nSnapshots = opt$n_snapshots,
    snapshotPhase = opt$snapshot_phase,
    asynchronousUpdates = isTRUE(opt$asynchronous_updates)
  )
  saveRDS(result, result_path)

  feature_loadings <- as.matrix(getFeatureLoadings(result))
  sample_factors <- as.matrix(getSampleFactors(result))
  feature_loadings <- ensure_pattern_matrix_orientation(feature_loadings, nrow(sce), gene_names)
  sample_factors <- ensure_pattern_matrix_orientation(sample_factors, ncol(sce), cell_names)

  colnames(feature_loadings) <- canonicalize_pattern_names(colnames(feature_loadings) %||% sprintf("Pattern%d", seq_len(ncol(feature_loadings))))
  colnames(sample_factors) <- canonicalize_pattern_names(colnames(sample_factors) %||% sprintf("Pattern%d", seq_len(ncol(sample_factors))))
  feature_loadings <- feature_loadings[, pattern_columns(colnames(feature_loadings)), drop = FALSE]
  sample_factors <- sample_factors[, colnames(feature_loadings), drop = FALSE]
  feature_loadings <- clean_matrix(feature_loadings)
  sample_factors <- clean_matrix(sample_factors)
  pattern_names <- colnames(feature_loadings)

  gene_loadings_csv <- file.path(outdir, sprintf("%s.gene_loadings.csv", stem))
  cell_scores_csv <- file.path(outdir, sprintf("%s.cell_scores.csv", stem))
  top_genes_csv <- file.path(outdir, sprintf("%s.top_genes.csv", stem))
  cell_metadata_patterns_csv <- file.path(outdir, sprintf("%s.discovery_cells_with_patterns.csv", stem))
  pattern_summary_csv <- file.path(outdir, sprintf("%s.pattern_summary.csv", stem))
  rq1_csv <- file.path(outdir, "RQ1_time_varying_patterns.csv")
  rq2_csv <- file.path(outdir, "RQ2_identity_vs_activity_patterns.csv")
  rq3_csv <- file.path(outdir, "RQ3_interferon_candidate_patterns.csv")
  means_by_time_csv <- file.path(outdir, "experimental_pattern_means_by_time.csv")
  means_by_time_celltype_csv <- file.path(outdir, "experimental_pattern_means_by_time_and_celltype.csv")

  write.csv(feature_loadings, gene_loadings_csv, quote = TRUE)
  write.csv(sample_factors, cell_scores_csv, quote = TRUE)
  top_gene_df <- top_genes_table(feature_loadings, opt$top_genes)
  write.csv(top_gene_df, top_genes_csv, row.names = FALSE, quote = TRUE)

  cell_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  cell_meta <- cell_meta[rownames(sample_factors), , drop = FALSE]
  cell_meta_patterns <- cbind(cell_meta, as.data.frame(sample_factors))
  write.csv(cell_meta_patterns, cell_metadata_patterns_csv, quote = TRUE)

  pattern_summary <- summarize_patterns(cell_meta_patterns, top_gene_df)
  write.csv(pattern_summary, pattern_summary_csv, row.names = FALSE, quote = TRUE)
  write.csv(pattern_summary[order(pattern_summary$eta_timepoint, decreasing = TRUE), ], rq1_csv, row.names = FALSE, quote = TRUE)
  rq2 <- pattern_summary[, c("pattern", "eta_timepoint", "eta_broad_cell_type", "pattern_class", "peak_timepoint")]
  rq2 <- rq2[order(rq2$pattern_class, -rq2$eta_timepoint), ]
  write.csv(rq2, rq2_csv, row.names = FALSE, quote = TRUE)
  write.csv(pattern_summary[pattern_summary$candidate_ifn_pattern, ], rq3_csv, row.names = FALSE, quote = TRUE)

  means_by_time <- pattern_means(cell_meta_patterns, pattern_names, c("timepoint_merged", "day_merged_numeric"))
  means_by_time <- means_by_time[order(means_by_time$day_merged_numeric), ]
  write.csv(means_by_time, means_by_time_csv, row.names = FALSE, quote = TRUE)
  means_by_time_cell <- pattern_means(cell_meta_patterns, pattern_names, c("broad_cell_type", "timepoint_merged", "day_merged_numeric"))
  means_by_time_cell <- means_by_time_cell[order(means_by_time_cell$broad_cell_type, means_by_time_cell$day_merged_numeric), ]
  write.csv(means_by_time_cell, means_by_time_celltype_csv, row.names = FALSE, quote = TRUE)

  md <- slot(result, "metadata")
  trace_df <- trace_table(md)
  write.csv(trace_df, trace_path, row.names = FALSE, quote = TRUE)
  diagnostics <- list(
    meanChiSq = md$meanChiSq %||% NA_real_,
    totalRunningTime = md$totalRunningTime %||% NA_real_,
    totalUpdates = md$totalUpdates %||% NA_real_,
    chisqHistory_length = length(md$chisq %||% numeric()),
    atomhistoryA_length = length(md$atomsA %||% numeric()),
    atomhistoryP_length = length(md$atomsP %||% numeric()),
    equilibrationSnapshotsA_length = length(md$equilibrationSnapshotsA %||% list()),
    equilibrationSnapshotsP_length = length(md$equilibrationSnapshotsP %||% list()),
    samplingSnapshotsA_length = length(md$samplingSnapshotsA %||% list()),
    samplingSnapshotsP_length = length(md$samplingSnapshotsP %||% list()),
    pumpStat_shape = safe_dim(md$pumpStat %||% NULL),
    pumpMatrix_shape = safe_dim(md$pumpMatrix %||% NULL),
    meanPatternAssignment_shape = safe_dim(md$meanPatternAssignment %||% NULL)
  )
  write_json_file(diagnostics_path, diagnostics)

  effective <- effective_pattern_summary(feature_loadings, sample_factors)
  redundancy <- redundancy_summary(feature_loadings, top_gene_df)
  pattern_counts <- table(pattern_summary$pattern_class)
  eta_series <- pattern_summary$eta_timepoint[is.finite(pattern_summary$eta_timepoint)]
  ifn_overlap_series <- pattern_summary$ifn_top_gene_overlap_top15[is.finite(pattern_summary$ifn_top_gene_overlap_top15)]
  ifn_rho_series <- pattern_summary$spearman_ifn_score_rho[is.finite(pattern_summary$spearman_ifn_score_rho)]
  top_genes_by_pattern <- lapply(pattern_names, function(pattern) top_gene_df$gene[top_gene_df$pattern == pattern])
  names(top_genes_by_pattern) <- pattern_names

  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  metrics <- c(metrics, diagnostics, redundancy, list(
    openmpSupport = openmp_support,
    runtime_sec = runtime_sec,
    completed_at_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
    orientation_msg = "R CoGAPS input used genes as rows and cells as columns from zellkonverter readH5AD.",
    gene_loadings_csv = gene_loadings_csv,
    cell_scores_csv = cell_scores_csv,
    top_genes_csv = top_genes_csv,
    cell_metadata_patterns_csv = cell_metadata_patterns_csv,
    pattern_summary_csv = pattern_summary_csv,
    rq1_time_varying_patterns_csv = rq1_csv,
    rq2_identity_vs_activity_patterns_csv = rq2_csv,
    rq3_interferon_candidate_patterns_csv = rq3_csv,
    experimental_pattern_means_by_time_csv = means_by_time_csv,
    experimental_pattern_means_by_time_and_celltype_csv = means_by_time_celltype_csv,
    pattern_names = pattern_names,
    n_patterns_effective = effective$n_effective,
    degenerate_patterns = effective$degenerate,
    effective_pattern_details = effective$rows,
    top_genes_by_pattern = top_genes_by_pattern,
    activity_like_pattern_count = table_count(pattern_counts, "activity-like"),
    identity_like_pattern_count = table_count(pattern_counts, "identity-like"),
    mixed_pattern_count = table_count(pattern_counts, "mixed"),
    unclassified_pattern_count = table_count(pattern_counts, "unclassified"),
    candidate_ifn_pattern_count = sum(pattern_summary$candidate_ifn_pattern, na.rm = TRUE),
    max_ifn_top_gene_overlap_top15 = if (length(ifn_overlap_series)) max(ifn_overlap_series) else 0,
    max_ifn_score_rho = if (length(ifn_rho_series)) max(ifn_rho_series) else 0,
    mean_eta_timepoint = if (length(eta_series)) mean(eta_series) else 0,
    max_eta_timepoint = if (length(eta_series)) max(eta_series) else 0
  ))
  write_json_file(metrics_path, metrics)
  message(sprintf("[DONE] %s finished in %.2f min", tag, runtime_sec / 60))
  message(sprintf("[DONE] metrics=%s", metrics_path))
}, error = function(e) {
  metrics$status <<- "error"
  metrics$runtime_sec <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  metrics$completed_at_utc <<- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  metrics$error_type <<- class(e)[[1]]
  metrics$error_message <<- conditionMessage(e)
  write_json_file(metrics_path, metrics)
  stop(e)
})
