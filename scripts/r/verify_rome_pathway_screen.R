# Independent, base-R check of held-out Rome results. Inputs are the source
# GSE25628 and GPL571 gzip files, plus a TSV projecting the public nominations
# onto the locked Entrez memberships (family, id, direction, gene; one row/gene).
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
read_gz_lines <- function(path) {
    stream <- gzfile(path, "rt")
    on.exit(close(stream))
    readLines(stream, warn = TRUE)
}
get_block <- function(lines, start, finish) {
    lo <- which(lines == start)
    hi <- which(lines == finish)
    stopifnot(length(lo) == 1L, length(hi) == 1L, hi - lo >= 2L)
    read.delim(text = lines[(lo + 1L):(hi - 1L)], check.names = FALSE,
               colClasses = "character", na.strings = character())
}
matrix_lines <- read_gz_lines(args[1L])
fields <- function(key) {
    entries <- matrix_lines[startsWith(matrix_lines, paste0(key, "\t"))]
    lapply(entries, function(line) scan(text = line, what = "", sep = "\t", quiet = TRUE)[-1L])
}
one_field <- function(key) {
    entries <- fields(key)
    stopifnot(length(entries) == 1L)
    entries[[1L]]
}
stopifnot(identical(one_field("!Series_geo_accession"), "GSE25628"),
          identical(one_field("!Series_platform_id"), "GPL571"))
clinical <- fields("!Sample_characteristics_ch1")
field_for <- function(prefix) {
    found <- Filter(function(x) all(startsWith(x, prefix)), clinical)
    stopifnot(length(found) == 1L)
    substring(found[[1L]], nchar(prefix) + 1L)
}
stopifnot(all(field_for("tissue: ") == "endometrium"))
status <- field_for("donor: ")
biopsy <- field_for("disease state: ")
case <- status == "patient" & biopsy == "Normal (eutopic)"
control <- status == "healthy donor" & biopsy == "Normal (control)"
lesion <- status == "patient" & biopsy == "Pathological (ectopic)"
stopifnot(sum(case) == 8, sum(control) == 6, sum(lesion) == 8, all(case | control | lesion))
eligible <- case | control
group <- case[eligible]

tab <- get_block(matrix_lines, "!series_matrix_table_begin", "!series_matrix_table_end")
stopifnot(identical(names(tab)[-1L], one_field("!Sample_geo_accession")),
          !anyDuplicated(tab[[1L]]), nrow(tab) == 22277L)
annotation <- get_block(read_gz_lines(args[2L]), "!platform_table_begin", "!platform_table_end")
stopifnot(!anyDuplicated(annotation$ID), setequal(annotation$ID, tab[[1L]]))
a <- annotation[match(tab[[1L]], annotation$ID), ]
valid <- !startsWith(a$ID, "AFFX") & grepl("^[1-9][0-9]*$", a[["Gene ID"]]) &
         grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", a[["Gene symbol"]])
genes <- a[["Gene ID"]]
stopifnot(sum(valid) == 19921L, length(unique(genes[valid])) == 12502L)
x <- do.call(cbind, lapply(tab[-1L][eligible], function(column) {
    tokens <- tolower(trimws(column))
    vals <- suppressWarnings(as.numeric(column))
    recognized <- tokens %in% c("", "null", "na", "nan", "+nan", "-nan")
    stopifnot(!any(is.na(vals) & !recognized))
    vals
}))
finite <- rowSums(!is.finite(x)) == 0L
stopifnot(all(finite[valid]), sum(finite) >= 20000L)
ranked <- apply(x[finite, , drop = FALSE], 2L, function(column) {
    100 * (rank(column, ties.method = "average") - 0.5) / length(column)
})
probe_ids <- genes[valid]
summed <- rowsum(ranked[valid[finite], , drop = FALSE], probe_ids)
genome <- summed / as.numeric(table(probe_ids)[rownames(summed)])

lock <- read.delim(args[3L], colClasses = "character", check.names = FALSE)
stopifnot(identical(names(lock), c("family", "id", "direction", "gene")),
          nrow(lock) > 0L, !anyDuplicated(lock[c("id", "gene")]),
          all(lock$family %in% c("gene", "pathway")), all(lock$gene %in% rownames(genome)),
          all(lock$direction %in% c("higher", "lower")),
          !any(lock$gene == "2354"), length(unique(lock$id)) <= 10L)
assignments <- combn(seq_along(group), sum(group))
contrast <- function(score, affected) mean(score[affected]) - mean(score[!affected])
out <- lapply(unique(lock$id), function(id) {
    rows <- lock[lock$id == id, ]
    stopifnot(length(unique(rows$family)) == 1L, length(unique(rows$direction)) == 1L)
    if (rows$family[1L] == "gene") {
        stopifnot(nrow(rows) == 1L, identical(id, rows$gene))
    } else {
        stopifnot(nrow(rows) >= 15L, nrow(rows) <= 200L)
    }
    y <- colMeans(genome[rows$gene, , drop = FALSE])
    delta <- contrast(y, group)
    sd_null <- sqrt(var(y) * (1 / sum(group) + 1 / sum(!group)))
    stopifnot(is.finite(sd_null))
    if (sd_null > 0) {
        tails <- if (rows$direction[1L] == "higher") 1 else -1
        distribution <- apply(assignments, 2L, function(ix) {
            (mean(y[ix]) - mean(y[-ix])) * tails / sd_null
        })
        extremes <- sum(distribution >= tails * delta / sd_null - 1e-10)
    } else {
        extremes <- ncol(assignments)
    }
    deletes <- vapply(seq_along(y), function(i) contrast(y[-i], group[-i]), numeric(1L))
    data.frame(id = id, difference = delta, extreme = extremes,
               one_sided_p = extremes / ncol(assignments),
               omit_low = min(deletes), omit_high = max(deletes))
})
result <- do.call(rbind, out)
result$holm <- p.adjust(result$one_sided_p, method = "holm")
stopifnot(ncol(assignments) == 3003L)
write.table(result, row.names = FALSE, quote = FALSE, sep = "\t")
