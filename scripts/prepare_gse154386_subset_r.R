#!/usr/bin/env Rscript

# Download, preprocess, and subset GSE154386 for Case Study 5.
# This optional full-reproduction script prepares large local objects, but it
# does not run CoGAPS.

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(jsonlite)
})

RAW_TAR_URL <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE154386&format=file"
RANDOM_SEED <- 42L
MIN_GENES <- 200L
MIN_COUNTS <- 500L
MAX_PCT_MT <- 20
MIN_CELLS_PER_GENE <- 10L
N_HVGS <- 5000L
DISCOVERY_MAX_CELLS_PER_SAMPLE <- 600L

BROAD_MARKERS <- list(
  Monocyte = c("LST1", "FCN1", "CTSS", "SAT1", "TYMP", "S100A8", "S100A9"),
  T_cell = c("CD3D", "CD3E", "TRAC", "LTB", "IL7R", "MALAT1"),
  NK_cell = c("NKG7", "GNLY", "PRF1", "CTSW", "KLRD1", "TYROBP"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "BANK1"),
  Plasmablast = c("MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"),
  Dendritic = c("FCER1A", "CST3", "CLEC10A", "CD1C", "HLA-DRA"),
  Neutrophil = c("FCGR3B", "CXCR2", "CSF3R", "MNDA", "S100A8", "S100A9")
)

PROGRAMS <- list(
  ifn_program = c(
    "IFITM1", "IFI6", "ISG15", "MX1", "IFIT1", "IFIT3",
    "IFI44L", "ISG20", "LY6E", "TRIM22", "OAS1", "OASL"
  ),
  plasmablast_program = c("MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"),
  translation_program = c("RPL4", "RPL5", "RPL6", "RPS3", "RPS8", "EEF2"),
  mito_program = c("MT-CYB", "MT-ND4", "MT-ND1", "MT-CO1", "MT-CO2")
)

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    workdir = "data/external/GSE154386",
    max_cells_per_sample = DISCOVERY_MAX_CELLS_PER_SAMPLE,
    n_hvgs = N_HVGS,
    seed = RANDOM_SEED
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- args[[i + 1L]]
    if (key == "--workdir") out$workdir <- value
    if (key == "--max-cells-per-sample") out$max_cells_per_sample <- as.integer(value)
    if (key == "--n-hvgs") out$n_hvgs <- as.integer(value)
    if (key == "--seed") out$seed <- as.integer(value)
    i <- i + 2L
  }
  out
}

download_if_needed <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message("[skip] ", dest, " already exists")
    return(invisible(dest))
  }
  message("[download] ", url, " -> ", dest)
  download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  invisible(dest)
}

extract_if_needed <- function(tar_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(out_dir, all.files = FALSE, recursive = FALSE)) > 0) {
    message("[skip] ", out_dir, " already has extracted files")
    return(invisible(out_dir))
  }
  message("[extract] ", tar_path, " -> ", out_dir)
  utils::untar(tar_path, exdir = out_dir)
  invisible(out_dir)
}

parse_sample_key <- function(sample_key) {
  match <- regexec("^(GSM[0-9]+)_([^_]+)_(D(?:neg)?[0-9]+)$", sample_key)
  parts <- regmatches(sample_key, match)[[1]]
  if (length(parts) == 0) stop("Could not parse sample key: ", sample_key)

  gsm_id <- parts[[2]]
  subject <- parts[[3]]
  tp <- parts[[4]]
  cohort <- if (startsWith(subject, "Subject")) "experimental" else "natural"
  day <- if (startsWith(tp, "Dneg")) -as.integer(sub("Dneg", "", tp)) else as.integer(sub("D", "", tp))

  data.frame(
    gsm_id = gsm_id,
    sample_key = sample_key,
    sample_id = gsm_id,
    subject = subject,
    cohort = cohort,
    timepoint_raw = tp,
    timepoint_display = paste0("D", day),
    day_numeric = day,
    is_reference = (cohort == "experimental" && day == 0) || (cohort == "natural" && day == 180),
    stringsAsFactors = FALSE
  )
}

strip_known_suffix <- function(fname) {
  suffix_map <- c(
    "_matrix.mtx.gz" = "matrix",
    "_matrix.mtx" = "matrix",
    "_barcodes.tsv.gz" = "barcodes",
    "_barcodes.tsv" = "barcodes",
    "_features.tsv.gz" = "features",
    "_features.tsv" = "features",
    "_genes.tsv.gz" = "genes",
    "_genes.tsv" = "genes"
  )
  for (suffix in names(suffix_map)) {
    if (endsWith(fname, suffix)) {
      return(list(key = substr(fname, 1L, nchar(fname) - nchar(suffix)), role = suffix_map[[suffix]]))
    }
  }
  list(key = fname, role = NA_character_)
}

discover_10x_triplets <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  groups <- list()
  for (path in files) {
    parsed <- strip_known_suffix(basename(path))
    if (is.na(parsed$role)) next
    groups[[parsed$key]][[parsed$role]] <- path
  }
  Filter(function(x) {
    !is.null(x$matrix) && !is.null(x$barcodes) && (!is.null(x$features) || !is.null(x$genes))
  }, groups)
}

read_table_gz <- function(path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, header = FALSE, stringsAsFactors = FALSE)
}

read_mtx_gz <- function(path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "rb") else file(path, "rb")
  on.exit(close(con), add = TRUE)
  Matrix::readMM(con)
}

make_unique <- function(x) make.unique(as.character(x), sep = "_")

read_sample_sce <- function(sample_key, triplet) {
  meta <- parse_sample_key(sample_key)
  mat <- read_mtx_gz(triplet$matrix) # features x cells
  barcodes <- read_table_gz(triplet$barcodes)[[1]]
  features_path <- triplet$features %||% triplet$genes
  features <- read_table_gz(features_path)

  gene_ids <- as.character(features[[1]])
  gene_symbols <- if (ncol(features) >= 2L) as.character(features[[2]]) else gene_ids
  feature_types <- if (ncol(features) >= 3L) as.character(features[[3]]) else rep("Gene Expression", length(gene_symbols))
  keep <- feature_types == "Gene Expression"
  if (!any(keep)) keep <- rep(TRUE, length(gene_symbols))

  mat <- as(mat[keep, , drop = FALSE], "dgCMatrix")
  gene_ids <- gene_ids[keep]
  gene_symbols <- make_unique(gene_symbols[keep])
  rownames(mat) <- gene_symbols
  colnames(mat) <- paste0(meta$sample_id, ":", barcodes)

  col_data <- meta[rep(1L, ncol(mat)), , drop = FALSE]
  col_data$barcode <- barcodes
  rownames(col_data) <- colnames(mat)

  sce <- SingleCellExperiment(
    assays = list(counts = mat),
    rowData = DataFrame(gene_id = gene_ids, gene_symbol = gene_symbols),
    colData = DataFrame(col_data)
  )
  sce
}

add_merged_timepoints <- function(sce) {
  day <- as.numeric(sce$day_numeric)
  label <- as.character(sce$timepoint_display)
  is_exp_d14_15 <- sce$cohort == "experimental" & day %in% c(14, 15)
  label[is_exp_d14_15] <- "D14/15"
  day[is_exp_d14_15] <- 14.5
  sce$timepoint_merged <- label
  sce$day_merged_numeric <- day
  sce
}

upper_outliers_by_sample <- function(values, sample_id, q = 0.995) {
  flags <- rep(FALSE, length(values))
  for (sample in unique(sample_id)) {
    idx <- which(sample_id == sample)
    if (length(idx) >= 10L) {
      threshold <- as.numeric(stats::quantile(values[idx], probs = q, na.rm = TRUE))
      flags[idx] <- values[idx] > threshold
    }
  }
  flags
}

score_gene_sets <- function(sce, gene_sets) {
  logmat <- logcounts(sce)
  for (label in names(gene_sets)) {
    present <- intersect(gene_sets[[label]], rownames(sce))
    if (length(present) >= 2L) {
      colData(sce)[[paste0(label, "_score")]] <- Matrix::colMeans(logmat[present, , drop = FALSE])
    }
  }
  sce
}

assign_broad_cell_type <- function(sce) {
  score_cols <- paste0(names(BROAD_MARKERS), "_score")
  score_cols <- score_cols[score_cols %in% colnames(colData(sce))]
  scores <- as.data.frame(colData(sce)[, score_cols, drop = FALSE])
  labels <- sub("_score$", "", score_cols)
  best <- labels[max.col(as.matrix(scores), ties.method = "first")]
  sce$broad_cell_type <- best
  sce
}

select_hvgs <- function(sce, n_hvgs) {
  logmat <- as.matrix(logcounts(sce))
  vars <- apply(logmat, 1L, stats::var)
  keep <- order(vars, decreasing = TRUE)[seq_len(min(n_hvgs, length(vars)))]
  sce[keep, ]
}

balanced_subset <- function(sce, max_cells_per_sample, seed) {
  if (is.na(max_cells_per_sample)) return(sce)
  set.seed(seed)
  keep <- unlist(lapply(split(seq_len(ncol(sce)), sce$sample_id), function(idx) {
    if (length(idx) <= max_cells_per_sample) idx else sample(idx, max_cells_per_sample)
  }), use.names = FALSE)
  sce[, keep]
}

write_h5ad_if_available <- function(sce, path) {
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    zellkonverter::writeH5AD(sce, path)
  } else {
    message("[skip] zellkonverter not installed, did not write ", path)
  }
}

main <- function() {
  args <- parse_args()
  workdir <- args$workdir
  raw_tar <- file.path(workdir, "GSE154386_RAW.tar")
  extract_dir <- file.path(workdir, "extracted")
  results_dir <- file.path(workdir, "results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  download_if_needed(RAW_TAR_URL, raw_tar)
  extract_if_needed(raw_tar, extract_dir)

  triplets <- discover_10x_triplets(extract_dir)
  if (length(triplets) == 0L) stop("No 10x matrix triplets found under ", extract_dir)

  sample_sces <- Map(read_sample_sce, names(triplets), triplets)
  sce_all <- do.call(cbind, sample_sces)
  sce_all <- add_merged_timepoints(sce_all)

  is_mito <- grepl("^MT-", rownames(sce_all), ignore.case = TRUE)
  qc <- scuttle::perCellQCMetrics(sce_all, subsets = list(Mito = which(is_mito)))
  sce_all$n_genes_by_counts <- qc$detected
  sce_all$total_counts <- qc$sum
  sce_all$pct_counts_mt <- qc$subsets_Mito_percent
  sce_all$high_genes_outlier <- upper_outliers_by_sample(qc$detected, sce_all$sample_id)
  sce_all$high_counts_outlier <- upper_outliers_by_sample(qc$sum, sce_all$sample_id)

  keep_cells <- sce_all$n_genes_by_counts >= MIN_GENES &
    sce_all$total_counts >= MIN_COUNTS &
    sce_all$pct_counts_mt <= MAX_PCT_MT &
    !sce_all$high_genes_outlier &
    !sce_all$high_counts_outlier
  sce_qc <- sce_all[, keep_cells]
  keep_genes <- Matrix::rowSums(counts(sce_qc) > 0) >= MIN_CELLS_PER_GENE
  sce_qc <- sce_qc[keep_genes, ]

  sce_qc <- scuttle::logNormCounts(sce_qc)
  sce_qc <- add_merged_timepoints(sce_qc)
  sce_qc <- score_gene_sets(sce_qc, BROAD_MARKERS)
  sce_qc <- score_gene_sets(sce_qc, PROGRAMS)

  sce_hvg <- select_hvgs(sce_qc, args$n_hvgs)
  sce_hvg <- assign_broad_cell_type(sce_hvg)
  sce_hvg <- add_merged_timepoints(sce_hvg)

  comp <- as.data.frame(colData(sce_hvg))
  comp <- aggregate(
    x = list(n_cells = rep(1L, nrow(comp))),
    by = comp[c("cohort", "timepoint_merged", "day_merged_numeric", "broad_cell_type")],
    FUN = sum
  )
  comp$fraction <- ave(comp$n_cells, comp$cohort, comp$timepoint_merged, FUN = function(x) x / sum(x))
  write.csv(comp, file.path(results_dir, "cell_composition_by_timepoint.csv"), row.names = FALSE)

  sce_exp <- sce_hvg[, sce_hvg$cohort == "experimental"]
  sce_nat <- sce_hvg[, sce_hvg$cohort == "natural"]
  sce_discovery <- balanced_subset(sce_exp, args$max_cells_per_sample, args$seed)
  sce_discovery <- add_merged_timepoints(sce_discovery)

  discovery_counts <- as.data.frame(colData(sce_discovery))
  discovery_counts <- aggregate(
    x = list(n_cells = rep(1L, nrow(discovery_counts))),
    by = discovery_counts[c("sample_id", "subject", "timepoint_merged", "day_merged_numeric", "broad_cell_type")],
    FUN = sum
  )
  write.csv(discovery_counts, file.path(results_dir, "experimental_discovery_cell_counts.csv"), row.names = FALSE)

  saveRDS(sce_hvg, file.path(workdir, "gse154386_preprocessed_hvg.rds"))
  saveRDS(sce_discovery, file.path(workdir, "gse154386_experimental_discovery_genes_x_cells.rds"))
  saveRDS(sce_nat, file.path(workdir, "gse154386_natural_projection_target.rds"))
  write_h5ad_if_available(sce_hvg, file.path(workdir, "gse154386_preprocessed_hvg.h5ad"))
  write_h5ad_if_available(sce_discovery, file.path(workdir, "gse154386_experimental_discovery_genes_x_cells.h5ad"))
  write_h5ad_if_available(sce_nat, file.path(workdir, "gse154386_natural_projection_target.h5ad"))

  manifest <- list(
    source = "GSE154386",
    raw_tar_url = RAW_TAR_URL,
    n_samples = length(triplets),
    all_cells_after_qc = ncol(sce_qc),
    hvg_shape = c(cells = ncol(sce_hvg), genes = nrow(sce_hvg)),
    experimental_discovery_shape = c(cells = ncol(sce_discovery), genes = nrow(sce_discovery)),
    natural_projection_shape = c(cells = ncol(sce_nat), genes = nrow(sce_nat)),
    max_cells_per_sample = args$max_cells_per_sample,
    seed = args$seed
  )
  write_json(manifest, file.path(results_dir, "preprocessing_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
}

main()
