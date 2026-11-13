# Fully synthetic microarray model tests; no network, patient records, or discoveries.
source("scripts/r/analyze_geo_microarray.R")
if (!requireNamespace("limma", quietly = TRUE)) stop("limma is required for this test")

root <- tempfile("synthetic_geo_")
dir.create(root)
set.seed(78)
sample_ids <- sprintf("GSM900%03d", 1:18)
phase <- rep(c("proliferative", "secretory"), each = 9)
group <- rep(rep(c("control", "case", "other"), each = 3), 2)
probe_ids <- sprintf("SYN_PROBE_%03d", 1:40)
mat <- matrix(rnorm(40 * 18, mean = 7, sd = 0.35), nrow = 40)
mat[1L, group == "case"] <- mat[1L, group == "case"] + 2
mat[2L, group == "other"] <- mat[2L, group == "other"] + 1.4
mat[40L, 1L] <- NA_real_
expression <- data.frame(ID_REF = probe_ids, mat, check.names = FALSE)
names(expression)[-1L] <- sample_ids
expression_path <- file.path(root, "expression.tsv")
utils::write.table(expression, expression_path, sep = "\t", row.names = FALSE,
                   quote = FALSE, na = "NA")
manifest <- data.frame(
  geo_accession = sample_ids, donor_id = paste0("SYN_DONOR_", seq_along(sample_ids)),
  include = "true", reviewed = "true", disease_group = group,
  phase = phase, tissue = "eutopic_endometrium",
  pathology = ifelse(group == "case", "endometriosis",
                     ifelse(group == "control", "none", "other_pelvic_condition")),
  hormone_treatment = "none", exclusion_reason = "",
  review_source_url = "https://example.org/synthetic-no-real-source",
  stringsAsFactors = FALSE
)
manifest_path <- file.path(root, "manifest.csv")
write_manifest <- function(x) utils::write.csv(x, manifest_path, row.names = FALSE, na = "")
write_manifest(manifest)
opts <- list(expression = expression_path, manifest = manifest_path,
             "out-dir" = file.path(root, "result"),
             "log2-confirmed" = "yes", "min-group" = 3L, seed = 20261112L)

run_microarray(opts)
case_result <- utils::read.csv(file.path(opts[["out-dir"]], "case_vs_control.csv"))
other_result <- utils::read.csv(file.path(opts[["out-dir"]], "other_vs_control.csv"))
null_result <- utils::read.csv(file.path(opts[["out-dir"]], "phase_permuted_control.csv"))
stopifnot(nrow(case_result) == 40L,
          case_result$logFC[1L] > 1.5,
          other_result$logFC[2L] > 0.8,
          case_result$status[40L] == "excluded_missing_intensity",
          is.na(case_result$P.Value[40L]),
          all(case_result$q_across_contrasts[1:39] >= case_result$P.Value[1:39]),
          nrow(null_result) == 40L,
          file.exists(file.path(opts[["out-dir"]], "sample_review_audit.csv")),
          any(grepl("n_included_independent_donors=18", readLines(file.path(opts[["out-dir"]], "analysis_report.txt")))))

expect_failure <- function(bad_manifest, pattern) {
  write_manifest(bad_manifest)
  rejected <- tryCatch({
    run_microarray(modifyList(opts, list("out-dir" = tempfile("blocked_geo_"))))
    FALSE
  }, error = function(e) grepl(pattern, conditionMessage(e), fixed = TRUE))
  if (!isTRUE(rejected)) stop("Expected guard failure: ", pattern)
}

bad <- manifest
bad$reviewed[1L] <- "false"
expect_failure(bad, "Every included sample must be reviewed")
bad <- manifest
bad$donor_id[2L] <- bad$donor_id[1L]
expect_failure(bad, "Non-independent samples")
bad <- manifest
bad$disease_group[1L] <- "healthy"
expect_failure(bad, "disease_group must be exactly")
bad <- manifest
bad$hormone_treatment[1L] <- "unknown"
expect_failure(bad, "hormone_treatment")
bad <- manifest
bad$hormone_treatment[1L] <- "hormonal_contraceptive"
expect_failure(bad, "separate protocol")
bad <- manifest
bad$phase[bad$disease_group == "case"] <- "case_only_phase"
expect_failure(bad, "no shared reviewed menstrual phase")

# A separate specificity contrast is conditional on reviewed other-pathology samples.
no_other <- manifest
no_other$include[no_other$disease_group == "other"] <- "false"
no_other$exclusion_reason[no_other$disease_group == "other"] <- "not_in_this_contrast"
write_manifest(no_other)
opts[["out-dir"]] <- file.path(root, "no_other")
run_microarray(opts)
stopifnot(!file.exists(file.path(opts[["out-dir"]], "other_vs_control.csv")),
          any(grepl("other_pathology_specificity=not_tested", readLines(file.path(opts[["out-dir"]], "analysis_report.txt")))))

bad_opts <- tryCatch({parse_args(c("--expression", "a", "--manifest", "b",
                                  "--out-dir", "c", "--log2-confirmed", "no")); FALSE},
                     error = function(e) grepl("log2", conditionMessage(e), fixed = TRUE))
stopifnot(bad_opts)
cat("Synthetic limma GEO microarray gates passed\n")
