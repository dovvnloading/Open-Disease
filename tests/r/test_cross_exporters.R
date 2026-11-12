# Both offline ingestion paths must agree for identical simulated raw counts.
# No Seurat RDS, human donors, genomic data or network required.
source("scripts/seurat_to_pseudobulk.R")
source("scripts/r/pseudobulk_core.R")

raw <- Matrix::sparseMatrix(
  i = c(1L, 2L, 3L, 2L, 1L, 3L),
  j = c(1L, 1L, 2L, 3L, 4L, 4L),
  x = c(5, 1, 2, 3, 4, 7),
  dims = c(3L, 4L),
  dimnames = list(c("SYN_G1", "SYN_G2", "SYN_G3"),
                  c("s1-c1", "s1-c2", "s2-c1", "s2-c2"))
)
meta <- data.frame(
  sample = c("sample1", "sample1", "sample2", "sample2"),
  AnnotationMain = rep("epithelial", 4L),
  row.names = colnames(raw), stringsAsFactors = FALSE
)
sample_map <- data.frame(
  sample_id = c("sample1", "sample2"),
  study_id = rep("toy_study", 2),
  donor_id = c("toy_case", "toy_control"),
  diagnosis = c("case", "control"),
  phase = rep("proliferative", 2),
  hormone_treatment = rep("none", 2),
  stringsAsFactors = FALSE
)
folder <- tempfile("open-disease-equivalence-")
dir.create(folder)
out_sample <- file.path(folder, "by_sample.csv")
out_audit <- file.path(folder, "by_sample_audit.csv")
out_cell <- file.path(folder, "by_cell.csv")
od_export_pseudobulk(raw, meta, sample_map, out_sample, out_audit,
                     min_cells = 1L)
cell_map <- cbind(cell_id = rownames(meta), sample_map[
  match(meta$sample, sample_map$sample_id),
  c("study_id", "donor_id", "diagnosis", "phase", "hormone_treatment"),
  drop = FALSE
], cell_type = meta$AnnotationMain)
cell_map <- as.data.frame(cell_map, stringsAsFactors = FALSE)
write_pseudobulk(raw, cell_map, out_cell)
a <- utils::read.csv(out_sample, stringsAsFactors = FALSE)
b <- utils::read.csv(out_cell, stringsAsFactors = FALSE)
order_cols <- c("study_id", "donor_id", "cell_type", "gene_id")
sort_rows <- function(frame) {
  frame[do.call(order, frame[order_cols]), names(frame), drop = FALSE]
}
a <- sort_rows(a)
b <- sort_rows(b)
rownames(a) <- NULL
rownames(b) <- NULL
stopifnot(isTRUE(all.equal(a, b, check.attributes = FALSE)),
          nrow(a) == 6L,
          all(a$n_cells == 2),
          all(a$library_size[a$donor_id == "toy_case"] == 8),
          all(a$library_size[a$donor_id == "toy_control"] == 14))
cat("PASS: both R exporters agree on all artificial donor/gene counts.\n")
