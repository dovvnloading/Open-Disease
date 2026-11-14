# Base R, with no Bioconductor packages required. Read only the named gene axis.
# The mapping script verifies the resulting label-order digest before use.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Expected epithelial RDS path and feature-label output path")
obj <- readRDS(args[[1L]])
ranges <- attr(obj, "rowRanges", exact = TRUE)
partitions <- attr(ranges, "partitioning", exact = TRUE)
labels <- attr(partitions, "NAMES", exact = TRUE)
stopifnot(is.character(labels), length(labels) == 29034L,
          !anyNA(labels), !anyDuplicated(labels), all(nzchar(labels)))
writeLines(labels, args[[2L]], useBytes = TRUE)
