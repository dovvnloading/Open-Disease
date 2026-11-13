#!/usr/bin/env Rscript
# Local-only analysis of reviewed GEO *log2 microarray* intensities using limma.
# Not a count model, patient classifier, gene-to-target map, or clinical result.

parse_args <- function(args) {
  if (length(args) %% 2L != 0L) stop("Pass --option value pairs", call. = FALSE)
  out <- list()
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || !is.null(out[[substring(key, 3L)]])) {
      stop("Unknown/duplicate option: ", key, call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
  }
  required <- c("expression", "manifest", "out-dir", "log2-confirmed")
  if (!all(required %in% names(out))) {
    stop("Required: --expression --manifest --out-dir --log2-confirmed yes", call. = FALSE)
  }
  if (!identical(out[["log2-confirmed"]], "yes")) {
    stop("Explicit log2 scale/source-processing review required", call. = FALSE)
  }
  unknown <- setdiff(names(out), c(required, "min-group", "seed"))
  if (length(unknown)) stop("Unknown options: ", paste(unknown, collapse = ", "), call. = FALSE)
  out[["min-group"]] <- if (is.null(out[["min-group"]])) 3L else as.integer(out[["min-group"]])
  out$seed <- if (is.null(out$seed)) 20261112L else as.integer(out$seed)
  if (is.na(out[["min-group"]]) || out[["min-group"]] < 3L || is.na(out$seed) || out$seed < 1L) {
    stop("min-group must be >= 3 and seed a positive integer", call. = FALSE)
  }
  out
}

parse_flag <- function(x, label) {
  flag <- tolower(trimws(as.character(x)))
  if (anyNA(flag) || !all(flag %in% c("true", "false"))) {
    stop(label, " must be explicitly true/false for every sample", call. = FALSE)
  }
  flag == "true"
}

require_text <- function(x, label) {
  if (anyNA(x) || any(!nzchar(trimws(as.character(x)))) ||
      any(tolower(trimws(as.character(x))) %in% c("unknown", "unreviewed", "n/a"))) {
    stop("Included samples lack independently reviewed ", label, call. = FALSE)
  }
}

make_table <- function(fit, contrast, probes, eligible, missing_n) {
  out <- data.frame(probe_id = probes,
                    status = ifelse(eligible, "tested", "excluded_missing_intensity"),
                    missing_samples = missing_n, stringsAsFactors = FALSE)
  cols <- c("logFC", "CI.L", "CI.R", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  for (col in cols) out[[col]] <- rep(NA_real_, length(probes))
  table <- limma::topTable(fit, coef = contrast, number = Inf,
                           sort.by = "none", confint = TRUE)
  if (!identical(rownames(table), probes[eligible])) {
    stop("limma output probe order does not match input", call. = FALSE)
  }
  for (col in intersect(cols, names(table))) out[[col]][eligible] <- table[[col]]
  out
}

run_microarray <- function(opts) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Install Bioconductor limma", call. = FALSE)
  manifest <- utils::read.csv(opts$manifest, stringsAsFactors = FALSE,
                              check.names = FALSE, na.strings = c("", "NA"))
  required <- c("geo_accession", "donor_id", "include", "reviewed", "disease_group",
                "phase", "tissue", "pathology", "hormone_treatment",
                "exclusion_reason", "review_source_url")
  if (!all(required %in% names(manifest))) {
    stop("Manifest missing columns: ", paste(setdiff(required, names(manifest)), collapse = ", "), call. = FALSE)
  }
  require_text(manifest$geo_accession, "GSM accession")
  if (anyDuplicated(manifest$geo_accession)) stop("Duplicate GEO sample accession", call. = FALSE)
  included <- parse_flag(manifest$include, "include")
  reviewed <- parse_flag(manifest$reviewed, "reviewed")
  if (any(included & !reviewed)) stop("Every included sample must be reviewed", call. = FALSE)
  if (any(!included & (is.na(manifest$exclusion_reason) |
                        !nzchar(trimws(manifest$exclusion_reason))))) {
    stop("Excluded samples require an explicit reason", call. = FALSE)
  }
  samples <- manifest[included, , drop = FALSE]
  if (!nrow(samples)) stop("Zero reviewed samples included", call. = FALSE)
  for (field in c("donor_id", "disease_group", "phase", "tissue", "pathology",
                  "hormone_treatment", "review_source_url")) {
    require_text(samples[[field]], field)
  }
  if (anyDuplicated(samples$donor_id)) {
    stop("Non-independent samples: repeated donor_id; paired model not implemented", call. = FALSE)
  }
  if (!all(samples$disease_group %in% c("control", "case", "other"))) {
    stop("disease_group must be exactly case/control/other", call. = FALSE)
  }
  if (!all(samples$tissue == "eutopic_endometrium")) {
    stop("Only reviewed eutopic_endometrium tissue is eligible", call. = FALSE)
  }
  if (!all(samples$hormone_treatment == "none")) {
    stop("Unknown or exposed hormone treatment requires a separate protocol", call. = FALSE)
  }
  if (!all(grepl("^https://", samples$review_source_url))) {
    stop("Included samples require an HTTPS clinical metadata review source", call. = FALSE)
  }
  counts <- table(samples$disease_group)
  min_group <- opts[["min-group"]]
  if (!all(c("case", "control") %in% names(counts)) ||
      any(counts[c("case", "control")] < min_group)) {
    stop("Insufficient independent case or control donors", call. = FALSE)
  }
  has_other <- "other" %in% names(counts)
  if (has_other && counts[["other"]] < min_group) {
    stop("Insufficient independent other-pathology donors", call. = FALSE)
  }
  phase_counts <- table(samples$phase, samples$disease_group)
  if (!any(phase_counts[, "control"] > 0 & phase_counts[, "case"] > 0)) {
    stop("Case/control groups have no shared reviewed menstrual phase", call. = FALSE)
  }
  if (has_other && !any(phase_counts[, "control"] > 0 & phase_counts[, "other"] > 0)) {
    stop("Other/control groups have no shared reviewed menstrual phase", call. = FALSE)
  }

  raw <- utils::read.delim(opts$expression, check.names = FALSE, stringsAsFactors = FALSE,
                           quote = '"', comment.char = "", na.strings = c("", "NA"))
  if (!"ID_REF" %in% names(raw) || ncol(raw) < 2L) stop("Missing ID_REF", call. = FALSE)
  probes <- as.character(raw$ID_REF)
  require_text(probes, "ID_REF")
  if (anyDuplicated(probes)) stop("Duplicate array probe ID", call. = FALSE)
  matrix_samples <- setdiff(names(raw), "ID_REF")
  if (anyDuplicated(names(raw)) || length(matrix_samples) != nrow(manifest) ||
      !setequal(matrix_samples, manifest$geo_accession)) {
    stop("Expression samples must match all manifest GSM accessions exactly", call. = FALSE)
  }
  intensity_frame <- raw[, samples$geo_accession, drop = FALSE]
  if (!all(vapply(intensity_frame, is.numeric, logical(1)))) {
    stop("All array intensities must be numeric or missing", call. = FALSE)
  }
  intensity <- as.matrix(intensity_frame)
  if (any(is.infinite(intensity)) || any(is.nan(intensity))) {
    stop("Non-finite array intensity detected", call. = FALSE)
  }
  missing_n <- rowSums(is.na(intensity))
  eligible <- missing_n == 0L  # conservative complete-case, outcome-blind probe gate
  if (sum(eligible) < 2L) stop("Fewer than two complete-case probes", call. = FALSE)
  mat <- intensity[eligible, , drop = FALSE]
  rownames(mat) <- probes[eligible]
  group_levels <- c("control", "case", if (has_other) "other")
  group <- factor(samples$disease_group, levels = group_levels)
  phase <- factor(samples$phase, levels = sort(unique(samples$phase)))
  formula <- if (nlevels(phase) > 1L) ~ 0 + group + phase else ~ 0 + group
  design <- stats::model.matrix(formula)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) {
    stop("Non-estimable group/phase model; review confounding and residual degrees of freedom", call. = FALSE)
  }
  contrast_names <- c("case_vs_control", if (has_other) "other_vs_control")
  contrasts <- matrix(0, nrow = ncol(design), ncol = length(contrast_names),
                      dimnames = list(colnames(design), contrast_names))
  contrasts["groupcase", "case_vs_control"] <- 1
  contrasts["groupcontrol", "case_vs_control"] <- -1
  if (has_other) {
    contrasts["groupother", "other_vs_control"] <- 1
    contrasts["groupcontrol", "other_vs_control"] <- -1
  }
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(mat, design), contrasts))
  results <- lapply(contrast_names, function(name) make_table(fit, name, probes, eligible, missing_n))
  names(results) <- contrast_names
  # Within-contrast q is limma's adj.P.Val; across planned contrasts is separate.
  all_p <- unlist(lapply(results, function(x) x$P.Value[eligible]), use.names = FALSE)
  all_q <- stats::p.adjust(all_p, method = "BH")
  k <- sum(eligible)
  for (i in seq_along(results)) {
    results[[i]]$q_across_contrasts <- rep(NA_real_, length(probes))
    results[[i]]$q_across_contrasts[eligible] <- all_q[((i - 1L) * k + 1L):(i * k)]
  }

  # Negative control: a fixed-seed, within-phase label permutation, NOT validation.
  set.seed(opts$seed)
  permuted <- as.character(group)
  for (level in levels(phase)) {
    index <- which(phase == level)
    before <- permuted[index]
    if (length(unique(before)) > 1L) {
      after <- sample(before, length(before), replace = FALSE)
      if (identical(after, before)) after <- c(before[-1L], before[1L])
      permuted[index] <- after
    }
  }
  if (identical(permuted, as.character(group))) {
    stop("Null-control permutation did not change any labels", call. = FALSE)
  }
  group <- factor(permuted, levels = group_levels)
  null_design <- stats::model.matrix(formula)
  if (!identical(colnames(null_design), colnames(design)) ||
      qr(null_design)$rank != ncol(null_design)) {
    stop("Negative-control design is not estimable", call. = FALSE)
  }
  null_fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(mat, null_design), contrasts))
  negative_control <- make_table(null_fit, "case_vs_control", probes, eligible, missing_n)

  out_dir <- opts[["out-dir"]]
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; refusing overwrite", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(results)) {
    utils::write.csv(results[[name]], file.path(out_dir, paste0(name, ".csv")),
                     row.names = FALSE, na = "NA")
  }
  utils::write.csv(negative_control, file.path(out_dir, "phase_permuted_control.csv"),
                   row.names = FALSE, na = "NA")
  utils::write.csv(manifest, file.path(out_dir, "sample_review_audit.csv"),
                   row.names = FALSE, na = "NA")
  utils::write.csv(as.data.frame(phase_counts), file.path(out_dir, "phase_group_counts.csv"),
                   row.names = FALSE)
  lines <- c(
    "scope=exploratory_GSE51981_probe_level_no_clinical_use",
    paste0("source_expression=", normalizePath(opts$expression)),
    paste0("source_manifest=", normalizePath(opts$manifest)),
    "assay=processed_log2_microarray_not_counts",
    paste0("n_manifest=", nrow(manifest)),
    paste0("n_included_independent_donors=", nrow(samples)),
    paste0("n_case=", counts[["case"]]),
    paste0("n_control=", counts[["control"]]),
    paste0("n_other=", if (has_other) counts[["other"]] else 0L),
    paste0("n_probes_total=", length(probes)),
    paste0("n_probes_tested_complete_case=", sum(eligible)),
    paste0("n_probes_excluded_missing=", sum(!eligible)),
    paste0("model=", paste(deparse(formula), collapse = "")),
    paste0("model_rank=", qr(design)$rank),
    paste0("residual_df=", nrow(design) - ncol(design)),
    paste0("contrasts=", paste(contrast_names, collapse = ",")),
    paste0("other_pathology_specificity=", if (has_other) "tested" else "not_tested_no_reviewed_other_samples"),
    "multiplicity=BH_within_each_contrast_and_across_all_planned_contrasts",
    paste0("negative_control=one_within_phase_label_permutation_seed_", opts$seed),
    paste0("negative_control_n_q_below_0.05=", sum(negative_control$adj.P.Val < 0.05, na.rm = TRUE)),
    paste0("limma_version=", as.character(utils::packageVersion("limma")))
  )
  writeLines(lines, file.path(out_dir, "analysis_report.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "session_info.txt"))
  invisible(lines)
}

if (sys.nframe() == 0L) {
  run_microarray(parse_args(commandArgs(trailingOnly = TRUE)))
}
