# Synthetic-only tests. Does not need SeuratObject, a human atlas or a network.
source("scripts/seurat_to_pseudobulk.R")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix must be available")

expect_failure <- function(expr, match_text) {
  encountered <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  stopifnot(is.character(encountered), grepl(match_text, encountered, fixed = TRUE))
}

toy <- Matrix::sparseMatrix(
  i = c(1L, 2L, 2L, 3L),
  j = c(1L, 2L, 3L, 4L),
  x = c(4, 3, 1, 6),
  dims = c(3L, 4L),
  dimnames = list(c("SYN1", "SYN2", "SYN3"),
                  c("CELL1", "CELL2", "CELL3", "CELL4")))
stopifnot(inherits(toy, "dgCMatrix"))
mapping <- data.frame(
  cell_id = c("CELL4", "CELL2", "CELL1", "CELL3"),
  study_id = rep("SYN_STUDY", 4),
  donor_id = c("control1", "case1", "case1", "control1"),
  diagnosis = c("control", "case", "case", "control"),
  phase = rep("proliferative", 4),
  hormone_treatment = rep("none", 4),
  cell_type = rep("epithelial", 4),
  stringsAsFactors = FALSE)
outfile <- tempfile(fileext = ".csv")
write_pseudobulk(toy, mapping, outfile)
output <- utils::read.csv(outfile, stringsAsFactors = FALSE)
stopifnot(nrow(output) == 6L, identical(sort(unique(output$donor_id)), c("case1", "control1")))
case <- output[output$donor_id == "case1", ]
control <- output[output$donor_id == "control1", ]
stopifnot(all(case$library_size == 7), all(control$library_size == 7))
stopifnot(all(case$n_cells == 2), all(control$n_cells == 2))
stopifnot(identical(as.numeric(case$counts[match(c("SYN1", "SYN2", "SYN3"), case$gene_id)]),
                    c(4, 3, 0)))
stopifnot(identical(as.numeric(control$counts[match(c("SYN1", "SYN2", "SYN3"), control$gene_id)]),
                    c(0, 1, 6)))
expect_failure(write_pseudobulk(toy, mapping, outfile), "Refusing to overwrite")
bad_counts <- toy
bad_counts@x[1L] <- 0.5
expect_failure(write_pseudobulk(bad_counts, mapping, tempfile(fileext = ".csv")),
               "nonnegative integer raw counts")
bad_map <- mapping
bad_map$phase[bad_map$cell_id == "CELL2"] <- "secretory"
expect_failure(write_pseudobulk(toy, bad_map, tempfile(fileext = ".csv")),
               "Inconsistent phase within donor")
bad_map <- mapping
bad_map$hormone_treatment[1L] <- "yes"
expect_failure(write_pseudobulk(toy, bad_map, tempfile(fileext = ".csv")),
               "documented absence of exogenous hormone")
bad_map <- mapping
bad_map$cell_id[1L] <- "CELL1"
expect_failure(write_pseudobulk(toy, bad_map, tempfile(fileext = ".csv")),
               "Duplicate metadata cell_id")
bad_map <- mapping[-1L, , drop = FALSE]
expect_failure(write_pseudobulk(toy, bad_map, tempfile(fileext = ".csv")),
               "exactly one row per RNA count cell")
unlink(outfile)
cat("PASS: sparse donor aggregation, explicit zeros, donor mapping, counts and consent-context guards (synthetic only)\n")
