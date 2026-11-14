# Independent aggregate cross-check of the frozen New Delhi rank contrast.
# Usage: Rscript scripts/r/verify_gse120103_fosb.R PATH_TO_PINNED_MATRIX
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
con <- gzfile(args[[1L]], "rt")
lines <- readLines(con, warn = TRUE)
close(con)
header_field <- function(key) {
    matching <- lines[startsWith(lines, paste0(key, "\t"))]
    stopifnot(length(matching) == 1L)
    scan(text = matching, what = "", sep = "\t", quiet = TRUE)[-1L]
}
begin <- which(lines == "!series_matrix_table_begin")
end <- which(lines == "!series_matrix_table_end")
stopifnot(length(begin) == 1L, length(end) == 1L, end > begin + 1L)
tab <- read.delim(text = lines[(begin + 1L):(end - 1L)], check.names = FALSE,
                  stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN", "null"))
stopifnot(identical(header_field("!Series_geo_accession"), "GSE120103"),
          identical(names(tab)[-1L], header_field("!Sample_geo_accession")),
          !anyDuplicated(names(tab)), !anyDuplicated(tab[[1L]]))
group <- sub(".*\\(Group ([12][AB])\\)$", "\\1", header_field("!Sample_source_name_ch1"))
stopifnot(identical(as.integer(table(group)), rep(9L, 4L)),
          setequal(unique(group), c("2A", "1A", "2B", "1B")))
values <- as.matrix(tab[-1L])
keep <- rowSums(!is.finite(values)) == 0
values <- values[keep, , drop = FALSE]
target <- which(tab[[1L]][keep] == "A_23_P429998")
stopifnot(length(target) == 1L, nrow(values) >= 20000L)
score <- apply(values, 2L, function(column) {
    100 * (rank(column, ties.method = "average")[target] - 0.5) / length(column)
})
estimate <- function(score, group) {
    means <- tapply(score, group, mean)
    c(fertile = unname(means["2A"] - means["1A"]),
      infertile = unname(means["2B"] - means["1B"]))
}
strata <- estimate(score, group)
effect <- mean(strata)
component <- tapply(score, group, function(x) var(x) / (4 * length(x)))
df <- sum(component)^2 / sum(component^2 / (as.numeric(table(group)) - 1L))
se <- sqrt(sum(component))
interval <- effect + qt(c(0.025, 0.975), df = df) * se
delete_one <- vapply(seq_along(score), function(i) mean(estimate(score[-i], group[-i])), numeric(1L))
print(c(total = nrow(tab), complete = sum(keep), strata, effect = effect, se = se,
        df = df, ci_low = interval[1L], ci_high = interval[2L],
        two_sided_p = 2 * pt(-abs(effect / se), df),
        loo_min = min(delete_one), loo_max = max(delete_one)), digits = 10)
