# Independent aggregate cross-check of the frozen RNA transcript contrast.
# Usage: Rscript scripts/r/verify_gse315857_fosb.R MATRIX_GZ SAMPLE_SOFT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
con <- gzfile(args[1L], "rt")
tab <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
close(con)
meta <- readLines(args[2L], warn = FALSE)
starts <- which(startsWith(meta, "^SAMPLE = "))
ends <- c(starts[-1L] - 1L, length(meta))
get_field <- function(block, tag) {
    lines <- block[startsWith(block, tag)]
    stopifnot(length(lines) == 1L)
    substring(lines, nchar(tag) + 1L)
}
samples <- lapply(seq_along(starts), function(i) {
    block <- meta[starts[i]:ends[i]]
    title <- get_field(block, "!Sample_title = ")
    stopifnot(grepl("^(Eosis|Control) Prolif", title))
    c(library = get_field(block, "!Sample_description = Library name: "),
      group = if (startsWith(title, "Eosis Prolif")) "case" else "control")
})
sample_data <- do.call(rbind, samples)
lib <- sample_data[, "library"]
stopifnot(length(lib) == 8L, !anyDuplicated(lib), !anyDuplicated(names(tab)),
          identical(sort(lib), sort(names(tab)[-1L])),
          identical(as.integer(table(sample_data[, "group"])), c(4L, 4L)))
reference <- tab[[1L]] != "" & !is.na(tab[[1L]])
ids <- tab[[1L]][reference]
stopifnot(sum(reference) == 22961L, sum(!reference) == 4L, !anyDuplicated(ids))
target <- which(ids == "ENST00000353609")
x <- as.matrix(tab[reference, lib])
stopifnot(length(target) == 1L, all(is.finite(x)), all(x >= 0))
pct <- apply(x, 2L, function(column) {
    score <- rank(column, ties.method = "average")[target]
    100 * (score - 1) / (length(column) - 1)
})
groups <- sample_data[, "group"]
medians <- function(x, g) tapply(x, g, median)
centers <- medians(pct, groups)
difference <- function(x, g) {
    m <- medians(x, g)
    unname(m["case"] - m["control"])
}
leave_one <- vapply(seq_along(pct), function(i) difference(pct[-i], groups[-i]), numeric(1L))
print(c(centers, difference = difference(pct, groups), loo_min = min(leave_one),
        loo_max = max(leave_one)), digits = 12)
