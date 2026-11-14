# Called by the Swiss analysis wrapper with checksum-pinned sources and a
# temporary private output. No printed expression values, cell or donor keys.
shared_glandular_donors <- function(meta, study, min_cells = 20L) {
    specimen <- as.character(meta$sample)
    subtype <- as.character(meta$AnnotationRefined)
    idx <- match(specimen, study$specimen)
    linked <- study[idx, , drop = FALSE]
    stopifnot(length(specimen) == length(subtype), !anyNA(subtype), !anyNA(idx),
              all(nzchar(specimen)), !anyNA(study$specimen), !anyDuplicated(study$specimen),
              all(linked$deg %in% c("True", "False")),
              all(linked$diagnosis %in% c("ENDO", "non-ENDO")),
              all(as.character(meta$EndometriosisStatus) == ifelse(linked$diagnosis == "ENDO", "Endo", "CTL")),
              all(as.character(meta$MenstrualCyclePhase) == linked$phase),
              all(linked$phase[linked$deg == "True"] == "proliferative"))
    glandular <- which(linked$deg == "True" & subtype == "glandular")
    sufficient <- which(tabulate(idx[glandular], nbins = nrow(study)) >= min_cells)
    capture <- study$capture[sufficient]
    case <- study$diagnosis[sufficient] == "ENDO"
    sequencing <- study$sequencing[sufficient]
    stopifnot(!anyNA(capture), all(nzchar(capture)), !anyNA(sequencing), all(nzchar(sequencing)),
              all(vapply(split(sequencing, capture), function(v) length(unique(v)) == 1L, logical(1))))
    shared <- intersect(capture[case], capture[!case])
    kept <- sufficient[capture %in% shared]
    groups <- lapply(kept, function(k) glandular[idx[glandular] == k])
    list(cells = groups, capture = match(study$capture[kept], sort(shared)),
         case = study$diagnosis[kept] == "ENDO", eligible = length(sufficient))
}

pathway_members <- function(mapping, features) {
    lock <- c("R-HSA-9008059", "R-HSA-5654741", "R-HSA-163560", "R-HSA-5654710", "R-HSA-432722")
    stopifnot(identical(names(mapping), c("family", "id", "entrez_id", "feature_label", "row_index")))
    path <- mapping[mapping$family == "pathway", , drop = FALSE]
    stopifnot(setequal(path$id, lock), all(grepl("^[1-9][0-9]*$", path$entrez_id)),
              !anyDuplicated(path[c("id", "entrez_id")]),
              all(nzchar(path$feature_label) == nzchar(path$row_index)),
              identical(tabulate(match(path$id, lock)), c(19L, 38L, 21L, 18L, 43L)))
    mapped <- path[nzchar(path$row_index), , drop = FALSE]
    rows <- as.integer(mapped$row_index)
    stopifnot(all(grepl("^[1-9][0-9]*$", mapped$row_index)), all(rows >= 1L & rows <= length(features)),
              identical(features[rows], mapped$feature_label), !anyDuplicated(mapped[c("id", "row_index")]),
              identical(tabulate(match(mapped$id, lock)), c(18L, 33L, 18L, 15L, 43L)))
    unique_mapping <- unique(mapped[c("entrez_id", "row_index")])
    stopifnot(!anyDuplicated(unique_mapping$entrez_id), !anyDuplicated(unique_mapping$row_index))
    split(rows, factor(mapped$id, levels = lock))
}

within_donor_rank <- function(totals) {
    # Ranking sums and means gives identical ties and order within one donor.
    100 * (rank(totals, ties.method = "average") - 0.5) / length(totals)
}

score_donor <- function(columns, n_features, pointer, index, value, members) {
    rows_needed <- sort(unique(unlist(members, use.names = FALSE)))
    totals <- numeric(n_features)
    accumulated <- numeric(length(rows_needed))
    for (col in columns) {
        rows <- counts <- integer(0)
        if (pointer[col] < pointer[col + 1L]) {
            pos <- seq.int(pointer[col] + 1L, pointer[col + 1L])
            rows <- index[pos] + 1L
            counts <- value[pos]
            stopifnot(!anyNA(rows), !is.unsorted(rows, strictly = TRUE), all(rows >= 1L & rows <= n_features),
                      all(is.finite(counts)), all(counts >= 0), all(counts == trunc(counts)))
            totals[rows] <- totals[rows] + counts
        }
        # For a value v, (midrank - 0.5) = (# rows below v) + (# rows equal v)/2.
        # Rank only stored positive entries, adding the exact implicit/explicit
        # zero block; requested missing entries receive the zero tie midrank.
        positive <- counts > 0
        positive_rows <- rows[positive]
        n_zero <- n_features - length(positive_rows)
        ranks <- rep(n_zero / 2, length(rows_needed))
        at <- match(rows_needed, positive_rows)
        seen <- !is.na(at)
        if (any(seen)) {
            nonzero_ranks <- rank(counts[positive], ties.method = "average")
            ranks[seen] <- n_zero + nonzero_ranks[at[seen]] - 0.5
        }
        accumulated <- accumulated + ranks
    }
    primary_ranks <- within_donor_rank(totals)
    cell_ranks <- accumulated * 100 / (n_features * as.double(length(columns)))
    result <- data.frame(endpoint = names(members),
                         primary = vapply(members, function(r) mean(primary_ranks[r]), numeric(1)),
                         cell_first = vapply(members, function(r) mean(cell_ranks[match(r, rows_needed)]), numeric(1)),
                         row.names = NULL)
    rbind(result, data.frame(endpoint = "GLOBAL", primary = mean(result$primary), cell_first = mean(result$cell_first)))
}

score_swiss_glandular <- function(rds, study_file, mapping_file, output) {
    object <- readRDS(rds)
    field <- function(x, name) attr(x, name, exact = TRUE)
    meta <- field(field(object, "colData"), "listData")
    cells <- field(field(object, "colData"), "rownames")
    features <- field(field(field(object, "rowRanges"), "partitioning"), "NAMES")
    assays <- field(field(field(object, "assays"), "data"), "listData")
    counts <- assays[["counts"]]
    stopifnot(identical(as.character(field(object, "class")), "cell_data_set"),
              setequal(names(assays), c("counts", "logcounts")),
              identical(as.character(field(counts, "class")), "dgCMatrix"),
              identical(field(counts, "Dim"), c(29034L, 20549L)),
              identical(field(counts, "Dimnames"), list(features, cells)),
              length(features) == 29034L, length(cells) == 20549L,
              !anyNA(features), all(nzchar(features)), !anyDuplicated(features),
              !anyNA(cells), all(nzchar(cells)), !anyDuplicated(cells),
              all(c("sample", "AnnotationRefined", "EndometriosisStatus", "MenstrualCyclePhase") %in% names(meta)),
              all(vapply(meta, length, integer(1)) == length(cells)))
    pointer <- field(counts, "p")
    index <- field(counts, "i")
    value <- field(counts, "x")
    stopifnot(length(pointer) == length(cells) + 1L, pointer[1L] == 0L,
              !is.unsorted(pointer), pointer[length(pointer)] == length(index),
              length(index) == length(value))
    read_private <- function(path) utils::read.delim(path, colClasses = "character", check.names = FALSE,
                                                     na.strings = character(), quote = "")
    study <- read_private(study_file)
    stopifnot(nrow(study) == 67L, length(unique(meta$sample)) == 58L,
              identical(names(study), c("specimen", "diagnosis", "phase", "capture", "sequencing", "deg")))
    donor <- shared_glandular_donors(meta, study)
    allocation <- table(donor$capture, donor$case)
    stopifnot(donor$eligible == 21L, length(donor$cells) == 13L, sum(donor$case) == 6L,
              nrow(allocation) == 4L, all(allocation > 0L),
              prod(choose(rowSums(allocation), allocation[, "TRUE"])) == 144L)
    members <- pathway_members(read_private(mapping_file), features)
    scores <- lapply(seq_along(donor$cells), function(k) {
        z <- score_donor(donor$cells[[k]], length(features), pointer, index, value, members)
        data.frame(donor_index = k, capture = donor$capture[k], case = as.integer(donor$case[k]), z)
    })
    out <- do.call(rbind, scores)
    stopifnot(nrow(out) == 78L, all(is.finite(out$primary)), all(is.finite(out$cell_first)))
    utils::write.table(out, output, quote = FALSE, row.names = FALSE, sep = "\t")
    cat(sprintf("glandular_cells=%d\n", sum(lengths(donor$cells))))
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    stopifnot(length(args) == 4L)
    do.call(score_swiss_glandular, as.list(args))
}
