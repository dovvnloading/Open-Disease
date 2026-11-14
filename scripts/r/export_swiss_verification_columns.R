# Copy the selected native CSC columns into a temporary private stream. All
# ranking, pseudo-bulk construction and randomization for this check use Python.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
object <- readRDS(args[1L])
slot_value <- function(x, name) attr(x, name, exact = TRUE)
col_data <- slot_value(object, "colData")
meta <- slot_value(col_data, "listData")
rows <- slot_value(slot_value(slot_value(object, "rowRanges"), "partitioning"), "NAMES")
assays <- slot_value(slot_value(slot_value(object, "assays"), "data"), "listData")
counts <- assays[["counts"]]
features <- length(rows)
cells <- slot_value(col_data, "rownames")
stopifnot(identical(slot_value(counts, "Dim"), c(29034L, 20549L)),
          identical(slot_value(counts, "Dimnames"), list(rows, cells)),
          !anyDuplicated(rows), !anyDuplicated(cells))
clin <- utils::read.delim(args[2L], colClasses = "character", check.names = FALSE)
clin_ix <- match(as.character(meta$sample), clin$specimen)
stopifnot(!anyNA(clin_ix), !anyDuplicated(clin$specimen),
          identical(as.character(meta$MenstrualCyclePhase), clin$phase[clin_ix]),
          all((meta$EndometriosisStatus == "Endo") == (clin$diagnosis[clin_ix] == "ENDO")))
rows_lock <- utils::read.delim(args[3L], colClasses = "character", na.strings = "")
mapped <- rows_lock$family == "pathway" & !is.na(rows_lock$row_index)
stopifnot(identical(rows[as.integer(rows_lock$row_index[mapped])], rows_lock$feature_label[mapped]))

# Donors must contribute at least 20 strictly glandular cells in this object
# and satisfy the authors' proliferative differential-expression inclusion rule.
select <- which(meta$AnnotationRefined == "glandular" & clin$deg[clin_ix] == "True" &
                clin$phase[clin_ix] == "proliferative")
size <- tabulate(clin_ix[select], nbins = nrow(clin))
enough <- which(size >= 20L)
by_group <- split(clin$capture[enough], clin$diagnosis[enough])
shared <- intersect(by_group[["ENDO"]], by_group[["non-ENDO"]])
donors <- enough[clin$capture[enough] %in% shared]
col_use <- select[clin_ix[select] %in% donors]
labels <- as.integer(clin$diagnosis[donors] == "ENDO")
batch <- as.integer(factor(clin$capture[donors]))
stopifnot(length(donors) == 13L, sum(labels) == 6L, length(unique(batch)) == 4L,
          all(tapply(clin$sequencing[donors], batch, function(x) length(unique(x)) == 1L)))
row_indices <- slot_value(counts, "i")
column_ptr <- slot_value(counts, "p")
values <- slot_value(counts, "x")
export <- function() {
    io <- file(args[4L], "wb")
    on.exit(close(io))
    put <- function(x) writeBin(as.integer(x), io, size = 4L, endian = "little")
    put(c(features, length(donors), length(col_use)))
    put(labels)
    put(batch)
    put(size[donors])
    donor_number <- match(clin_ix, donors) - 1L
    for (j in col_use) {
        count <- column_ptr[j + 1L] - column_ptr[j]
        put(c(donor_number[j], count))
        if (count == 0L) next
        positions <- column_ptr[j] + seq_len(count)
        r <- row_indices[positions]
        v <- values[positions]
        stopifnot(all(r >= 0L & r < features), !is.unsorted(r, strictly = TRUE),
                  all(is.finite(v)), all(v >= 0 & v == trunc(v) & v <= .Machine$integer.max))
        put(r)
        put(v)
    }
}
export()
