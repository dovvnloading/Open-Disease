# Pure Matrix tests: do not require a real Seurat object or human data.
source("scripts/r/pseudobulk_core.R")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package missing")
mat <- Matrix::Matrix(
  matrix(c(5, 0, 1, 0, 3, 2, 2, 1, 0), nrow = 3,
         dimnames = list(c("GENE_A", "GENE_B", "GENE_C"),
                         c("cell_a", "cell_b", "cell_c"))),
  sparse = TRUE
)
metadata <- data.frame(
  sample = c("s1", "s1", "s2"),
  AnnotationMain = c("epithelial", "epithelial", "epithelial"),
  row.names = c("cell_a", "cell_b", "cell_c"),
  stringsAsFactors = FALSE
)
mapping <- data.frame(
  sample_id = c("s1", "s2"), study_id = c("synthetic", "synthetic"),
  donor_id = c("D1", "D2"), diagnosis = c("case", "control"),
  phase = c("proliferative", "proliferative"),
  hormone_treatment = c("none", "none"),
  stringsAsFactors = FALSE
)
dir <- tempfile("od-r-")
dir.create(dir)
out <- file.path(dir, "pseudobulk.csv")
audit <- file.path(dir, "audit.csv")
on.exit <- NULL
result <- od_export_pseudobulk(mat, metadata, mapping, out, audit, min_cells = 1L)
pb <- utils::read.csv(out, stringsAsFactors = FALSE)
qc <- utils::read.csv(audit, stringsAsFactors = FALSE)
stopifnot(
  result$groups == 2L, result$gene_count == 3L,
  nrow(pb) == 6L, nrow(qc) == 2L,
  all(pb$library_size[pb$donor_id == "D1"] == 11),
  all(pb$n_cells[pb$donor_id == "D1"] == 2L),
  pb$counts[pb$donor_id == "D1" & pb$gene_id == "GENE_B"] == 3,
  pb$counts[pb$donor_id == "D2" & pb$gene_id == "GENE_B"] == 1
)
# Refuse to overwrite original controlled results.
err_overwrite <- tryCatch({
  od_export_pseudobulk(mat, metadata, mapping, out, audit, min_cells = 1L)
  FALSE
}, error = function(e) TRUE)
stopifnot(err_overwrite)
# Reject a broken sample-to-donor mapping; fail before writing.
bad_map <- mapping
bad_map$diagnosis[[1]] <- "positive"
err_labels <- tryCatch({
  od_export_pseudobulk(mat, metadata, bad_map,
                       file.path(dir, "bad.csv"), file.path(dir, "bad_audit.csv"),
                       min_cells = 1L)
  FALSE
}, error = function(e) TRUE)
stopifnot(err_labels)
# Missing cells in metadata may never silently disappear.
bad_meta <- metadata[-1, , drop = FALSE]
err_cell <- tryCatch({
  od_export_pseudobulk(mat, bad_meta, mapping,
                       file.path(dir, "bad2.csv"), file.path(dir, "bad2_audit.csv"),
                       min_cells = 1L)
  FALSE
}, error = function(e) TRUE)
stopifnot(err_cell)
cat("R pseudobulk adapter tests passed (100% artificial counts).\n")
