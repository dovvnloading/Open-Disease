# Offline, structure-only inspection of the checksum-pinned GEO epithelial RDS.
# Run via scripts/audit_gse266265_epithelial.py with an isolated base-R process.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
object <- readRDS(args[1L])
field <- function(x, key) attr(x, key, exact = TRUE)
stopifnot(identical(as.character(field(object, "class")), "cell_data_set"))
meta <- field(field(object, "colData"), "listData")
cells <- field(field(object, "colData"), "rownames")
ids <- field(field(field(object, "rowRanges"), "partitioning"), "NAMES")
assays <- field(field(field(object, "assays"), "data"), "listData")
alt <- field(field(field(object, "int_colData"), "listData")[["altExps"]], "listData")
stopifnot(setequal(names(assays), c("counts", "logcounts")), !length(alt),
          !anyNA(cells), all(nzchar(cells)), !anyDuplicated(cells),
          !anyNA(ids), all(nzchar(ids)), !anyDuplicated(ids),
          all(vapply(meta, length, integer(1)) == length(cells)))

# Examine storage integrity and numeric type over bounded chunks without
# extracting or reporting any individual cell, participant, gene or abundance.
for (name in names(assays)) {
    m <- assays[[name]]
    shape <- field(m, "Dim")
    stopifnot(identical(as.character(field(m, "class")), "dgCMatrix"),
              identical(shape, c(length(ids), length(cells))),
              identical(field(m, "Dimnames"), list(ids, cells)))
    colptr <- field(m, "p")
    index <- field(m, "i")
    values <- field(m, "x")
    stopifnot(length(colptr) == shape[2] + 1L, colptr[1L] == 0L,
              !any(diff(colptr) < 0), colptr[length(colptr)] == length(index),
              length(index) == length(values))
    all_integer <- TRUE
    all_nonnegative <- TRUE
    for (start in seq.int(1L, length(values), by = 100000L)) {
        batch <- seq.int(start, min(length(values), start + 99999L))
        v <- values[batch]
        stopifnot(all(is.finite(v)), all(index[batch] >= 0L & index[batch] < shape[1L]))
        all_integer <- all_integer && all(v == trunc(v))
        all_nonnegative <- all_nonnegative && all(v >= 0)
    }
    cat(sprintf("assay=%s class=dgCMatrix genes=%d cells=%d finite=TRUE nonnegative=%s integer_valued=%s\n",
                name, shape[1L], shape[2L], all_nonnegative, all_integer))
}
cat("row_labels_matching_ENSG=", sum(grepl("^ENSG[0-9]+([.][0-9]+)?$", ids)), "\n", sep = "")
cat("row_labels_symbol_format=", sum(grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", ids)), "\n", sep = "")
targets <- c("ABO", "KIF3C", "ZNF835", "GLTP", "F11")
cat("published_candidate_symbols_present=", paste(targets[targets %in% ids], collapse = ","), "\n", sep = "")
cat("published_candidate_symbols_without_exact_row_label=", paste(setdiff(targets, ids), collapse = ","), "\n", sep = "")
cat("alternate_experiments=", length(alt), "\n", sep = "")
cat("metadata_columns=", paste(names(meta), collapse = ","), "\n", sep = "")
stopifnot(all(c("sample", "EndometriosisStatus", "EndometriosisGrade",
                "MenstrualCyclePhase", "MenstrualCyclePhase_main",
                "AnnotationMain", "AnnotationRefined") %in% names(meta)))
specimen <- as.character(meta$sample)
stopifnot(!anyNA(specimen), all(nzchar(specimen)))
clinical <- data.frame(status = meta$EndometriosisStatus, grade = meta$EndometriosisGrade,
                       phase = as.character(meta$MenstrualCyclePhase),
                       broad_phase = as.character(meta$MenstrualCyclePhase_main))
stopifnot(!anyNA(clinical), all(clinical$status %in% c("CTL", "Endo")),
          all(clinical$grade %in% c("CTL", "mild", "severe")),
          all(clinical$phase %in% c("periovulatory", "proliferative", "secretory_early",
                                    "secretory_mid", "secretory_late")),
          all(clinical$broad_phase %in% c("periovulatory", "proliferative", "secretory")),
          all(meta$AnnotationMain[!is.na(meta$AnnotationMain)] == "epithelial"),
          all(meta$AnnotationRefined[!is.na(meta$AnnotationRefined)] %in% c("glandular", "glandular early-secretory",
                "lumenal", "lumenal 1", "lumenal 2", "mid-secretory", "TP63+/KRT5+")))
per_sample <- unique(data.frame(specimen, clinical))
if (anyDuplicated(per_sample$specimen))
    stop("Native sample label has incompatible disease, grade or phase records")
cat("unique_samples=", nrow(per_sample), "\n", sep = "")
cat("nonempty_refined_epithelial_subtypes=", length(unique(meta$AnnotationRefined[!is.na(meta$AnnotationRefined)])), "\n", sep = "")
cat("cells_missing_main_label=", sum(is.na(meta$AnnotationMain)), "\n", sep = "")
cat("cells_missing_refined_label=", sum(is.na(meta$AnnotationRefined)), "\n", sep = "")
cat("samples_without_refined_labels=", sum(tapply(!is.na(meta$AnnotationRefined), specimen, sum) == 0), "\n", sep = "")
cat("sample disease by phase (counts of unique sample labels; cell counts do not enter)\n")
print(with(per_sample, table(status, phase)))
cat("samples with at least 100 epithelial cells=", sum(table(specimen) >= 100L), "\n", sep = "")

study <- utils::read.delim(args[2L], colClasses = "character", na.strings = "", check.names = FALSE)
ix <- match(per_sample$specimen, study$specimen)
stopifnot(nrow(study) == 67L, !anyNA(study$specimen), !anyDuplicated(study$specimen),
          !anyNA(ix), !anyDuplicated(ix))
linked <- study[ix, ]
stopifnot(all(linked$diagnosis %in% c("ENDO", "non-ENDO")),
          all(ifelse(linked$diagnosis == "ENDO", "Endo", "CTL") == per_sample$status),
          all(linked$phase == per_sample$phase), !anyNA(linked$capture), !anyNA(linked$sequencing),
          all(linked$deg %in% c("True", "False")))
cat("source_linked_samples=", nrow(linked), "\n", sep = "")
cat("linked_capture_batches=", length(unique(linked$capture)), "\n", sep = "")
cat("linked_sequencing_batches=", length(unique(linked$sequencing)), "\n", sep = "")
strict <- linked[linked$deg == "True", , drop = FALSE]
stopifnot(all(strict$phase == "proliferative"))
cat("full_authors_DEG_population_eligible=")
print(table(study$diagnosis[study$deg == "True"]))
cat("linked_authors_DEG_population_eligible=")
print(table(strict$diagnosis))
for (scope in c("proliferative", "authors_DEG")) {
    subset <- if (scope == "authors_DEG") strict else linked[linked$phase == "proliferative", ]
    totals <- table(subset$capture, subset$diagnosis)
    both <- apply(totals > 0, 1L, all)
    cat(sprintf("%s_capture_common_support: batches=%d cases=%d controls=%d\n",
                scope, sum(both), sum(totals[both, "ENDO"]), sum(totals[both, "non-ENDO"])))
}
