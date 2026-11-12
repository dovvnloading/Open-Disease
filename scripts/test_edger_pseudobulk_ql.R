# Fake SYN### gene data only. No patient or disease measurements.
source("scripts/edger_pseudobulk_ql.R")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(error) || !grepl(pattern, error, fixed = TRUE)) {
    stop(paste("Expected error containing", pattern, "but got:", error), call. = FALSE)
  }
}

input <- "results/synthetic-demo/pseudobulk.csv"
if (!file.exists(input)) stop("Run bash scripts/run_demo.sh first")
s <- read_stratum(input, "synthetic-discovery", "proliferative", "epithelial")
stopifnot(s$n_case == 6, s$n_control == 6, nrow(s$counts) == 50,
          ncol(s$counts) == 12)
res <- run_ql(s, adjust = "none")
stopifnot(nrow(res) == 50, all(res$status %in% c("tested", "filtered_low_expression")))
marker <- res[res$gene_id == "SYN000", , drop = FALSE]
stopifnot(nrow(marker) == 1L, marker$status == "tested",
          marker$logFC_case_vs_control > 1, marker$q_value < 0.10,
          marker$q_value >= marker$p_value)

raw <- utils::read.csv(input, stringsAsFactors = FALSE)
bad <- raw
bad$library_size[1L] <- bad$library_size[1L] + 20L
path <- tempfile(fileext = ".csv")
utils::write.csv(bad, path, row.names = FALSE)
expect_error(read_stratum(path, "synthetic-discovery", "proliferative", "epithelial"),
             "Donor sample metadata conflict")
unlink(path)

bad <- raw
ix <- which(bad$study_id == "synthetic-discovery" &
            bad$phase == "proliferative" & bad$cell_type == "epithelial")[1L]
bad <- bad[-ix, , drop = FALSE]
path <- tempfile(fileext = ".csv")
utils::write.csv(bad, path, row.names = FALSE)
expect_error(read_stratum(path, "synthetic-discovery", "proliferative", "epithelial"),
             "Incomplete donor x gene matrix")
unlink(path)

bad <- raw
bad$counts[1L] <- 1.2
path <- tempfile(fileext = ".csv")
utils::write.csv(bad, path, row.names = FALSE)
expect_error(read_stratum(path, "synthetic-discovery", "proliferative", "epithelial"),
             "exact nonnegative integers")
unlink(path)

samples <- s$samples
cov <- data.frame(
  study_id = "synthetic-discovery", donor_id = samples$donor_id,
  batch = rep(c("A", "B"), length.out = nrow(samples)),
  age_years = rep(c(27, 35, 41, 48, 30, 38), length.out = nrow(samples)),
  stringsAsFactors = FALSE)
covfile <- tempfile(fileext = ".csv")
utils::write.csv(cov, covfile, row.names = FALSE)
adj <- run_ql(s, covariate_csv = covfile, adjust = "batch+age")
stopifnot(nrow(adj) == 50, all(adj$adjustment == "batch+age"),
          adj$logFC_case_vs_control[adj$gene_id == "SYN000"] > 1)

expect_error(run_ql(s, adjust = "batch+age"), "covariate CSV is required")
expect_error(run_ql(s, covariate_csv = covfile, adjust = "none"),
             "do not silently ignore")
wrong <- cov
wrong$batch <- ifelse(samples$diagnosis == "case", "case-only", "control-only")
wrongfile <- tempfile(fileext = ".csv")
utils::write.csv(wrong, wrongfile, row.names = FALSE)
expect_error(run_ql(s, covariate_csv = wrongfile, adjust = "batch+age"),
             "rank deficient")
expect_error(run_ql(s, min_residual_df = 20), "Insufficient residual degrees")
unlink(c(covfile, wrongfile))
cat("PASS: edgeR QL synthetic effects; donor, design, count and covariate failure gates\n")
