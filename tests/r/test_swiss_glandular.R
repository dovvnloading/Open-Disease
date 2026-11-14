# Base-R regression checks. All arrays, labels and keys below are synthetic.
source("scripts/r/score_swiss_glandular.R")

# Features include all-cell zeros, ties at nonzero levels, an empty cell, and
# an explicitly stored zero. The pathways overlap and have unequal lengths.
dense <- matrix(c(0, 0, 1, 3, 1, 0, 4,
                  0, 0, 0, 0, 0, 0, 0,
                  0, 1, 2, 0, 0, 0, 4,
                  0, 0, 0, 0, 0, 0, 0), nrow = 7)
rows <- list(c(1, 3, 4, 5, 7), integer(0), c(2, 3, 7), 6)
i <- as.integer(unlist(rows, use.names = FALSE) - 1L)
x <- unlist(lapply(seq_along(rows), function(k) dense[rows[[k]], k]), use.names = FALSE)
p <- c(0L, cumsum(lengths(rows)))
members <- list(A = c(1L, 3L), B = c(2L, 3L, 4L))
cols <- 1:4
result <- score_donor(cols, nrow(dense), p, i, x, members)
percent <- function(v) (rank(v, ties.method = "average") - 0.5) * 100 / length(v)
primary <- percent(rowMeans(dense))
cell <- rowMeans(apply(dense, 2L, percent))
expected <- function(values) {
    pathway <- vapply(members, function(ix) mean(values[ix]), numeric(1))
    unname(c(pathway, mean(pathway)))
}
stopifnot(identical(result$endpoint, c("A", "B", "GLOBAL")),
          isTRUE(all.equal(result$primary, expected(primary), tolerance = 1e-12)),
          isTRUE(all.equal(result$cell_first, expected(cell), tolerance = 1e-12)),
          any(abs(result$primary - result$cell_first) > 1),
          all(score_donor(2L, 7L, p, i, x, members)$cell_first == 50),
          inherits(try(score_donor(cols, 7L, p, replace(i, 2L, i[1L]), x, members), silent = TRUE), "try-error"))

# Clinical joins reference source labels explicitly, permit non-eligible rows
# with missing clinical data and require both groups in each included batch.
study <- data.frame(specimen = as.character(LETTERS[1:7]),
                    diagnosis = c("ENDO", "non-ENDO", "ENDO", "ENDO", "non-ENDO", "non-ENDO", "non-ENDO"),
                    phase = c(rep("proliferative", 5), "secretory_mid", "NA"),
                    capture = c("red", "red", "blue", "red", "red", "red", "red"),
                    sequencing = c("one", "one", "two", "one", "one", "one", "one"),
                    deg = c(rep("True", 5), "False", "NA"))
link <- c(1, 3, 2, 1, 3, 2, 4, 5, 5, 6, 6)
meta <- list(sample = study$specimen[link],
             AnnotationRefined = c(rep("glandular", 8), "lumenal", "glandular", "glandular"),
             EndometriosisStatus = ifelse(study$diagnosis[link] == "ENDO", "Endo", "CTL"),
             MenstrualCyclePhase = study$phase[link])
design <- shared_glandular_donors(meta, study, 2)
stopifnot(design$eligible == 3L, identical(design$cells, list(c(1L, 4L), c(3L, 6L))),
          identical(design$case, c(TRUE, FALSE)), identical(design$capture, c(1L, 1L)))
meta$EndometriosisStatus[1L] <- "CTL"
stopifnot(inherits(try(shared_glandular_donors(meta, study, 2), silent = TRUE), "try-error"))
cat("Swiss glandular synthetic checks passed\n")
