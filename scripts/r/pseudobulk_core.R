# Open-Disease independently authored Matrix/Seurat count aggregation.
# NEVER copy private participant data into the public repository.
#
# Input: full feature-by-cell raw RNA count matrix, cell metadata with
# row names matching cell IDs, and a separately audited local sample-to-donor
# mapping. Output: complete donor x cell type x gene pseudobulk CSV.
# This file uses public package APIs only, and contains no third-party
# unlicensed scientific implementation.

od_export_pseudobulk <- function(counts, cell_metadata, sample_map,
                                 output_csv, audit_csv,
                                 sample_col = "sample",
                                 celltype_col = "AnnotationMain",
                                 min_cells = 20L) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("R package Matrix is required")
  }
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    stop("counts must be a numeric gene-by-cell matrix or Matrix object")
  }
  if (!is.data.frame(cell_metadata) || !is.data.frame(sample_map)) {
    stop("cell_metadata and sample_map must be data.frames")
  }
  if (length(dim(counts)) != 2L || nrow(counts) < 1L || ncol(counts) < 1L) {
    stop("count matrix must contain genes and cells")
  }
  if (length(min_cells) != 1L || is.na(min_cells) ||
      min_cells != floor(min_cells) || min_cells < 1L) {
    stop("min_cells must be a positive integer")
  }
  if (length(output_csv) != 1L || length(audit_csv) != 1L ||
      normalizePath(dirname(output_csv), mustWork = FALSE) ==
      normalizePath(dirname(audit_csv), mustWork = FALSE) &&
      basename(output_csv) == basename(audit_csv)) {
    stop("Output and audit must be different paths")
  }
  if (file.exists(output_csv) || file.exists(audit_csv)) {
    stop("Refusing to overwrite an existing output or audit file")
  }
  gene_ids <- rownames(counts)
  cell_ids <- colnames(counts)
  if (is.null(gene_ids) || anyNA(gene_ids) ||
      any(!nzchar(trimws(gene_ids))) || anyDuplicated(gene_ids)) {
    stop("Unique, nonblank gene IDs are required")
  }
  if (is.null(cell_ids) || anyNA(cell_ids) ||
      any(!nzchar(trimws(cell_ids))) || anyDuplicated(cell_ids)) {
    stop("Unique, nonblank cell IDs are required")
  }
  if (is.null(rownames(cell_metadata)) ||
      anyDuplicated(rownames(cell_metadata)) ||
      !setequal(cell_ids, rownames(cell_metadata))) {
    stop("Cell metadata row names must match all count-matrix cell IDs exactly")
  }
  cell_metadata <- cell_metadata[cell_ids, , drop = FALSE]
  if (!all(c(sample_col, celltype_col) %in% names(cell_metadata))) {
    stop("Configured sample and cell-type metadata columns are missing")
  }
  required <- c("sample_id", "study_id", "donor_id", "diagnosis",
                "phase", "hormone_treatment")
  if (!all(required %in% names(sample_map))) {
    stop(paste("sample_map requires columns:", paste(required, collapse = ", ")))
  }
  for (field in required) {
    x <- sample_map[[field]]
    if (anyNA(x) || any(!nzchar(trimws(as.character(x))))) {
      stop(paste("Blank or NA in mapping column", field))
    }
  }
  if (anyDuplicated(as.character(sample_map$sample_id))) {
    stop("Sample IDs must be unique within the local mapping")
  }
  if (!all(sample_map$diagnosis %in% c("case", "control"))) {
    stop("diagnosis must be explicitly normalized to case/control")
  }
  if (!all(sample_map$phase %in% c("proliferative", "secretory"))) {
    stop("First baseline accepts only proliferative and secretory phase")
  }
  if (!all(sample_map$hormone_treatment == "none")) {
    stop("First baseline requires hormone_treatment=none; define treated cohort separately")
  }
  sample_id <- as.character(cell_metadata[[sample_col]])
  cell_type <- as.character(cell_metadata[[celltype_col]])
  if (anyNA(sample_id) || anyNA(cell_type) ||
      any(!nzchar(trimws(sample_id))) || any(!nzchar(trimws(cell_type)))) {
    stop("Missing sample or cell-type metadata")
  }
  m <- match(sample_id, as.character(sample_map$sample_id))
  if (anyNA(m)) {
    stop("Some cells have no audited sample-to-donor mapping")
  }
  joined <- sample_map[m, required, drop = FALSE]
  joined$cell_type <- cell_type
  rownames(joined) <- cell_ids

  # A participant may contribute more than one sample, but diagnosis,
  # cycle-phase, and treatment must not change within the (study, donor).
  donors <- unique(joined[, c("study_id", "donor_id"), drop = FALSE])
  for (i in seq_len(nrow(donors))) {
    selected <- joined$study_id == donors$study_id[i] &
                joined$donor_id == donors$donor_id[i]
    attributes <- joined[selected, c("diagnosis", "phase", "hormone_treatment"), drop = FALSE]
    if (any(vapply(attributes, function(x) length(unique(x)) > 1L, logical(1)))) {
      stop("Donor diagnosis, phase, or treatment conflicts across samples")
    }
  }

  # Inspect nonzero entries of sparse matrices without materializing the
  # enormous dense expression matrix. Require actual raw integer counts.
  values <- if (inherits(counts, "sparseMatrix")) {
    if (!("x" %in% slotNames(counts))) stop("Raw RNA count matrix must be numeric")
    methods::slot(counts, "x")
  } else {
    as.vector(counts)
  }
  if (!is.numeric(values) || anyNA(values) ||
      any(!is.finite(values)) || any(values < 0) ||
      any(values != floor(values))) {
    stop("Counts must contain nonnegative finite raw integers")
  }

  grouping <- unique(joined[, c("study_id", "donor_id",
                                  "diagnosis", "phase", "hormone_treatment", "cell_type"),
                            drop = FALSE])
  grouping <- grouping[do.call(order, grouping), , drop = FALSE]
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(audit_csv), recursive = TRUE, showWarnings = FALSE)

  # Both outputs are written to temporary files and renamed upon successful
  # completion, avoiding a misleading partially exported 'completed' file.
  tmp <- tempfile("pseudobulk-", tmpdir = dirname(output_csv), fileext = ".csv")
  tmp_audit <- tempfile("pseudobulk-audit-", tmpdir = dirname(audit_csv), fileext = ".csv")
  on.exit(unlink(c(tmp, tmp_audit)), add = TRUE)
  con <- file(tmp, open = "wt")
  con_open <- TRUE
  on.exit(if (con_open) close(con), add = TRUE)
  fields <- c("study_id", "donor_id", "diagnosis", "phase",
              "hormone_treatment", "cell_type", "gene_id",
              "counts", "library_size", "n_cells")
  writeLines(paste(fields, collapse = ","), con = con)

  n_written <- 0L
  auditing <- vector("list", nrow(grouping))
  for (i in seq_len(nrow(grouping))) {
    g <- grouping[i, , drop = FALSE]
    idx <- which(joined$study_id == g$study_id &
                 joined$donor_id == g$donor_id &
                 joined$cell_type == g$cell_type)
    n <- length(idx)
    status <- if (n < min_cells) "excluded_low_cell_count" else "exported"
    library_size <- NA_real_
    if (n >= min_cells) {
      aggregated <- Matrix::rowSums(counts[, idx, drop = FALSE])
      library_size <- sum(aggregated)
      if (!is.finite(library_size) || library_size <= 0) {
        stop("Cannot export a zero or invalid library size")
      }
      block <- data.frame(
        study_id = rep(g$study_id, length(gene_ids)),
        donor_id = rep(g$donor_id, length(gene_ids)),
        diagnosis = rep(g$diagnosis, length(gene_ids)),
        phase = rep(g$phase, length(gene_ids)),
        hormone_treatment = rep(g$hormone_treatment, length(gene_ids)),
        cell_type = rep(g$cell_type, length(gene_ids)),
        gene_id = gene_ids,
        counts = as.numeric(aggregated),
        library_size = rep(library_size, length(gene_ids)),
        n_cells = rep(n, length(gene_ids)),
        stringsAsFactors = FALSE
      )
      utils::write.table(block, con, append = TRUE, quote = TRUE,
                         sep = ",", col.names = FALSE, row.names = FALSE,
                         na = "")
      n_written <- n_written + 1L
    }
    auditing[[i]] <- data.frame(
      study_id = g$study_id,
      donor_id = g$donor_id,
      diagnosis = g$diagnosis,
      phase = g$phase,
      cell_type = g$cell_type,
      n_cells = n,
      library_size = library_size,
      status = status,
      stringsAsFactors = FALSE
    )
  }
  if (n_written == 0L) stop("No strata reached min_cells: nothing exported")

  close(con)
  con_open <- FALSE
  utils::write.csv(do.call(rbind, auditing), tmp_audit, row.names = FALSE)
  if (!file.rename(tmp, output_csv)) stop("Could not finalize pseudobulk CSV")
  if (!file.rename(tmp_audit, audit_csv)) {
    unlink(output_csv)
    stop("Could not finalize audit CSV")
  }
  message("Exported ", n_written, " donor/cell-type aggregates. ")
  message("If derived from humans, BOTH outputs are controlled research data, not safe for GitHub.")
  invisible(list(groups = n_written, gene_count = length(gene_ids),
                 excluded = nrow(grouping) - n_written))
}
