# Sparse Seurat RNA-count export to donor x cell-type x gene pseudobulk.
# NONCLINICAL RESEARCH. Inputs and outputs may be human-derived: keep local.
# Requires SeuratObject >=5 and Matrix for the live RDS adapter.
# The aggregation algorithm itself is tested with artificial sparse matrices.

stop_clean <- function(message) stop(message, call. = FALSE)

validate_raw_counts <- function(counts) {
  if (!inherits(counts, "dgCMatrix")) {
    stop_clean("Expected an unnormalized dgCMatrix RNA counts layer. Do not pass integrated/log-normalized data.")
  }
  if (nrow(counts) == 0L || ncol(counts) == 0L) stop_clean("Empty count matrix")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop_clean("Raw matrix must have gene row names and cell column names")
  }
  if (anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(trimws(rownames(counts)))) || any(!nzchar(trimws(colnames(counts))))) {
    stop_clean("Missing or blank gene/cell identifiers")
  }
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop_clean("Duplicate gene or cell identifiers in count matrix")
  }
  values <- counts@x
  if (anyNA(values) || any(!is.finite(values)) ||
      any(values < 0) || any(values != floor(values))) {
    stop_clean("RNA count matrix must contain finite nonnegative integer raw counts")
  }
  invisible(TRUE)
}

validate_export_metadata <- function(counts, metadata) {
  required <- c("cell_id", "study_id", "donor_id", "diagnosis",
                "phase", "hormone_treatment", "cell_type")
  if (!all(required %in% names(metadata))) {
    stop_clean(paste("Missing metadata columns:", paste(setdiff(required, names(metadata)), collapse = ", ")))
  }
  if (anyDuplicated(names(metadata))) stop_clean("Duplicate metadata column names")
  if (nrow(metadata) != ncol(counts)) {
    stop_clean("Metadata must contain exactly one row per RNA count cell; missing/extra rows")
  }
  for (name in required) {
    column <- as.character(metadata[[name]])
    if (anyNA(column) || any(!nzchar(trimws(column)))) {
      stop_clean(paste("Missing or blank metadata:", name))
    }
    # Prevent compound-key collisions, without guessing identifiers.
    if (any(grepl("\u001f", column, fixed = TRUE))) {
      stop_clean(paste("Prohibited control character in metadata:", name))
    }
    metadata[[name]] <- column
  }
  if (anyDuplicated(metadata$cell_id)) stop_clean("Duplicate metadata cell_id")
  positions <- match(colnames(counts), metadata$cell_id)
  if (anyNA(positions)) stop_clean("Metadata cell_id values do not match count matrix cell names")
  metadata <- metadata[positions, , drop = FALSE]
  if (!all(metadata$diagnosis %in% c("case", "control"))) {
    stop_clean("diagnosis must be manually mapped to case or control, not inferred")
  }
  if (!all(metadata$hormone_treatment == "none")) {
    stop_clean("Baseline only supports documented absence of exogenous hormone treatment")
  }
  donor_key <- interaction(metadata$study_id, metadata$donor_id,
                           drop = TRUE, sep = "\u001f")
  for (name in c("diagnosis", "phase", "hormone_treatment")) {
    inconsistent <- vapply(split(metadata[[name]], donor_key),
                           function(v) length(unique(v)) != 1L, logical(1))
    if (any(inconsistent)) stop_clean(paste("Inconsistent", name, "within donor"))
  }
  metadata
}

write_pseudobulk <- function(counts, metadata, out_path) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop_clean("Matrix package required")
  validate_raw_counts(counts)
  metadata <- validate_export_metadata(counts, metadata)
  if (file.exists(out_path)) stop_clean("Refusing to overwrite an existing output file")
  if (!dir.exists(dirname(out_path))) stop_clean("Output directory must already exist")

  # One donor/cell-type group at a time. Never densify the full atlas.
  keys <- interaction(metadata$study_id, metadata$donor_id, metadata$cell_type,
                      drop = TRUE, sep = "\u001f")
  groups <- split(seq_len(nrow(metadata)), keys)
  genes <- rownames(counts)
  tempfile_path <- tempfile(pattern = ".pseudobulk-", tmpdir = dirname(out_path), fileext = ".csv")
  on.exit(if (file.exists(tempfile_path)) unlink(tempfile_path), add = TRUE)
  group_number <- 0L
  for (idx in groups) {
    raw_counts <- as.numeric(Matrix::rowSums(counts[, idx, drop = FALSE]))
    library_size <- sum(raw_counts)
    if (!is.finite(library_size) || library_size <= 0) {
      stop_clean("Sample has a zero or invalid whole-transcriptome library size")
    }
    m <- metadata[idx[1L], , drop = FALSE]
    n_genes <- length(genes)
    frame <- data.frame(
      study_id = rep(m$study_id, n_genes),
      donor_id = rep(m$donor_id, n_genes),
      diagnosis = rep(m$diagnosis, n_genes),
      phase = rep(m$phase, n_genes),
      hormone_treatment = rep(m$hormone_treatment, n_genes),
      cell_type = rep(m$cell_type, n_genes),
      gene_id = genes,
      counts = raw_counts,
      library_size = rep(library_size, n_genes),
      n_cells = rep(length(idx), n_genes),
      stringsAsFactors = FALSE, check.names = FALSE)
    first <- group_number == 0L
    utils::write.table(frame, file = tempfile_path, sep = ",", row.names = FALSE,
                       col.names = first, append = !first, quote = TRUE, qmethod = "double")
    group_number <- group_number + 1L
  }
  if (!file.rename(tempfile_path, out_path)) {
    stop_clean("Could not atomically move pseudobulk output into place")
  }
  message("Wrote ", group_number, " donor/cell-type groups; ", length(genes),
          " gene rows each. Human-derived output: access-controlled storage only.")
  invisible(out_path)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  help <- paste(
    "NONCLINICAL local Seurat RNA-count to pseudobulk export.",
    "Requires a manually reviewed cell mapping CSV with:",
    "cell_id,study_id,donor_id,diagnosis,phase,hormone_treatment,cell_type",
    "Example:",
    "Rscript scripts/seurat_to_pseudobulk.R --rds /approved/GSE266265_EndoAtlas.rds --metadata /approved/cell_mapping.csv --out /approved/pseudobulk.csv",
    "Never publish individual-derived files or call synthetic results biological discoveries.",
    sep = "\n")
  if (length(args) == 1L && args[[1L]] == "--help") {
    cat(help, "\n")
    return(invisible(NULL))
  }
  if (length(args) != 6L || anyDuplicated(args[c(1L, 3L, 5L)]) ||
      !setequal(args[c(1L, 3L, 5L)], c("--rds", "--metadata", "--out"))) {
    stop_clean(paste("Exactly --rds PATH --metadata CSV --out CSV required.\n", help))
  }
  opts <- setNames(args[c(2L, 4L, 6L)], args[c(1L, 3L, 5L)])
  if (!file.exists(opts[["--rds"]]) || !file.exists(opts[["--metadata"]])) {
    stop_clean("Missing local input file")
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop_clean("Install SeuratObject in your approved R environment; not downloaded by this script")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) stop_clean("Matrix package required")
  message("READING RDS: use a high-memory, approved analysis environment")
  object <- readRDS(opts[["--rds"]])
  raw <- tryCatch(
    SeuratObject::GetAssayData(object, assay = "RNA", layer = "counts"),
    error = function(e) stop_clean(paste(
      "Could not obtain a unique raw RNA counts layer (no automatic layer joining):",
      conditionMessage(e)))
  )
  metadata <- utils::read.csv(opts[["--metadata"]], colClasses = "character",
                              check.names = FALSE, na.strings = c("NA", ""))
  write_pseudobulk(raw, metadata, opts[["--out"]])
}

if (sys.nframe() == 0L) main()
