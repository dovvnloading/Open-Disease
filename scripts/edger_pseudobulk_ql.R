# Open-Disease exploratory edgeR v4 donor-level pseudobulk quasi-likelihood analysis.
# NOT CLINICAL. No data acquisition or access is performed here.
stop_clean <- function(message) stop(message, call. = FALSE)

required <- function(frame, columns) {
  missing <- setdiff(columns, names(frame))
  if (length(missing)) stop_clean(paste("Missing input columns:", paste(missing, collapse = ", ")))
  if (anyDuplicated(names(frame))) stop_clean("Duplicate CSV column names")
}

integer_counts <- function(x, field, positive = FALSE) {
  n <- suppressWarnings(as.numeric(x))
  if (anyNA(n) || any(!is.finite(n)) || any(n < 0) ||
      any(n > 9007199254740991) || any(n != floor(n))) {
    stop_clean(paste(field, "must contain exact nonnegative integers"))
  }
  if (positive && any(n <= 0)) stop_clean(paste(field, "must be positive"))
  n
}

read_stratum <- function(path, study, phase, cell_type,
                         min_donors = 4L, min_cells = 20L, min_genes = 20L) {
  if (min_donors < 4L || min_cells < 1L || min_genes < 20L) {
    stop_clean("Minimum software gates: 4 donors/group, 1 cell/sample, 20 genes")
  }
  if (!file.exists(path)) stop_clean("Input pseudobulk CSV not found")
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  ids <- c("study_id", "donor_id", "diagnosis", "phase",
           "hormone_treatment", "cell_type", "gene_id")
  required(df, c(ids, "counts", "library_size", "n_cells"))
  if (nrow(df) == 0L) stop_clean("Empty pseudobulk CSV")
  for (column in ids) {
    x <- as.character(df[[column]])
    if (anyNA(x) || any(!nzchar(trimws(x)))) stop_clean(paste("Missing identifier:", column))
  }
  if (!all(df$diagnosis %in% c("case", "control"))) {
    stop_clean("Diagnosis must be explicitly case or control")
  }
  if (!all(df$hormone_treatment == "none")) {
    stop_clean("Untreated-only baseline: missing/treatment statuses cannot be inferred")
  }
  df$counts <- integer_counts(df$counts, "counts")
  df$library_size <- integer_counts(df$library_size, "library_size", TRUE)
  df$n_cells <- integer_counts(df$n_cells, "n_cells", TRUE)
  # A participant cannot change phase, diagnosis, or treatment elsewhere in
  # the full input, even in data excluded by the requested comparison.
  context <- unique(df[, c("study_id", "donor_id", "diagnosis",
                           "phase", "hormone_treatment")])
  if (anyDuplicated(context[, c("study_id", "donor_id")])) {
    stop_clean("Donor metadata conflict across dataset")
  }
  df <- df[df$study_id == study & df$phase == phase &
           df$cell_type == cell_type, , drop = FALSE]
  if (!nrow(df)) stop_clean("Study/phase/cell type not present")
  samples <- unique(df[, c("donor_id", "diagnosis", "library_size", "n_cells")])
  if (anyDuplicated(samples$donor_id)) stop_clean("Donor sample metadata conflict")
  if (anyDuplicated(df[, c("donor_id", "gene_id")])) {
    stop_clean("Duplicate donor/gene record")
  }
  donors <- sort(samples$donor_id)
  genes <- sort(unique(df$gene_id))
  if (length(genes) < min_genes) stop_clean("Too few assayed genes")
  if (nrow(df) != length(donors) * length(genes)) {
    stop_clean("Incomplete donor x gene matrix; missing observations cannot be imputed as zero")
  }
  samples <- samples[match(donors, samples$donor_id), , drop = FALSE]
  if (any(samples$n_cells < min_cells)) stop_clean("Too few cells per donor/cell-type sample")
  n_case <- sum(samples$diagnosis == "case")
  n_control <- sum(samples$diagnosis == "control")
  if (min(n_case, n_control) < min_donors) {
    stop_clean("Insufficient independent case/control donors")
  }
  counts <- matrix(0, nrow = length(genes), ncol = length(donors),
                   dimnames = list(genes, donors))
  counts[cbind(match(df$gene_id, genes), match(df$donor_id, donors))] <- df$counts
  # Stronger condition than the Python schematic validator: only full
  # raw transcriptome exports, not selected candidate-gene subsets.
  if (!isTRUE(all.equal(unname(colSums(counts)),
                        as.numeric(samples$library_size), tolerance = 0))) {
    stop_clean("library_size must equal the full exported assay gene count sum")
  }
  list(counts = counts, samples = samples, study = study, phase = phase,
       cell_type = cell_type, n_case = n_case, n_control = n_control)
}

donor_covariates <- function(samples, csv, study, adjust) {
  choices <- c("none", "batch", "age", "batch+age")
  if (!adjust %in% choices) stop_clean("adjust must be none, batch, age, or batch+age")
  if (adjust == "none") {
    if (!is.null(csv)) stop_clean("Covariates supplied but adjust=none: do not silently ignore them")
    return(samples)
  }
  if (is.null(csv) || !file.exists(csv)) {
    stop_clean("A donor-level covariate CSV is required for adjusted models")
  }
  df <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  fields <- c("study_id", "donor_id")
  if (grepl("batch", adjust, fixed = TRUE)) fields <- c(fields, "batch")
  if (grepl("age", adjust, fixed = TRUE)) fields <- c(fields, "age_years")
  required(df, fields)
  if (anyNA(df$study_id) || anyNA(df$donor_id) ||
      anyDuplicated(df[, c("study_id", "donor_id")])) {
    stop_clean("Invalid/duplicated donor keys in covariate CSV")
  }
  df <- df[df$study_id == study, , drop = FALSE]
  lookup <- match(samples$donor_id, df$donor_id)
  if (anyNA(lookup)) stop_clean("Missing matched donor covariates")
  if ("batch" %in% fields) {
    x <- as.character(df$batch[lookup])
    if (anyNA(x) || any(!nzchar(trimws(x)))) stop_clean("Missing batch labels")
    samples$batch <- factor(x)
    if (nlevels(samples$batch) < 2L) stop_clean("Single batch level: choose a prespecified alternate design")
  }
  if ("age_years" %in% fields) {
    x <- suppressWarnings(as.numeric(df$age_years[lookup]))
    if (anyNA(x) || any(!is.finite(x)) || any(x <= 0)) stop_clean("Invalid age_years values")
    samples$age_years <- x
  }
  samples
}

run_ql <- function(stratum, covariate_csv = NULL, adjust = "none",
                   min_residual_df = 4L, min_tested_genes = 20L) {
  if (min_residual_df < 4 || min_tested_genes < 20) {
    stop_clean("Minimum gates require residual df >=4 and >=20 tested genes")
  }
  if (!requireNamespace("edgeR", quietly = TRUE) ||
      packageVersion("edgeR") < "4.0.0") {
    stop_clean("A matching Bioconductor edgeR version >=4.0 is required")
  }
  samples <- donor_covariates(stratum$samples, covariate_csv,
                              stratum$study, adjust)
  samples$diagnosis <- factor(samples$diagnosis, levels = c("control", "case"))
  terms <- "diagnosis"
  if (grepl("batch", adjust, fixed = TRUE)) terms <- c(terms, "batch")
  if (grepl("age", adjust, fixed = TRUE)) terms <- c(terms, "age_years")
  design <- stats::model.matrix(stats::reformulate(terms), data = samples)
  if (!("diagnosiscase" %in% colnames(design))) {
    stop_clean("Could not construct case/control contrast")
  }
  if (qr(design)$rank < ncol(design)) {
    stop_clean("Non-estimable design: rank deficient or fully confounded")
  }
  residual_df <- nrow(design) - ncol(design)
  if (residual_df < min_residual_df) {
    stop_clean("Insufficient residual degrees of freedom")
  }
  if (!identical(colnames(stratum$counts), as.character(samples$donor_id))) {
    stop_clean("Counts and donor metadata order misaligned")
  }
  y <- edgeR::DGEList(counts = stratum$counts)
  keep <- edgeR::filterByExpr(y, design = design)
  if (sum(keep) < min_tested_genes) {
    stop_clean("Too few genes survive edgeR filterByExpr")
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  fit <- edgeR::glmQLFit(y, design = design, legacy = FALSE, robust = FALSE)
  test <- edgeR::glmQLFTest(fit, coef = which(colnames(design) == "diagnosiscase"))
  table <- test$table
  if (!all(c("logFC", "logCPM", "F", "PValue") %in% colnames(table))) {
    stop_clean("Unexpected edgeR test statistics")
  }
  if (anyNA(table$PValue) || any(!is.finite(table$PValue))) {
    stop_clean("Undefined edgeR p values")
  }
  all_genes <- rownames(stratum$counts)
  result <- data.frame(
    study_id = stratum$study, phase = stratum$phase, cell_type = stratum$cell_type,
    gene_id = all_genes, n_case = stratum$n_case, n_control = stratum$n_control,
    logFC_case_vs_control = NA_real_, logCPM = NA_real_, quasi_F = NA_real_,
    p_value = NA_real_, q_value = NA_real_,
    status = ifelse(keep, "tested", "filtered_low_expression"),
    adjustment = adjust, residual_df = residual_df,
    method = "edgeR_v4_QL_exploratory",
    edgeR_version = as.character(packageVersion("edgeR")),
    stringsAsFactors = FALSE)
  pos <- match(rownames(table), all_genes)
  if (anyNA(pos) || length(pos) != sum(keep)) stop_clean("edgeR gene IDs inconsistent")
  result$logFC_case_vs_control[pos] <- table$logFC
  result$logCPM[pos] <- table$logCPM
  result$quasi_F[pos] <- table$F
  result$p_value[pos] <- table$PValue
  result$q_value[pos] <- stats::p.adjust(table$PValue, method = "BH")
  result[order(result$gene_id), , drop = FALSE]
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  usage <- paste(
    "EXPLORATORY edgeR QL, NOT CLINICAL, NOT a validated endometriosis discovery.",
    "Required: --pseudobulk CSV --study STUDY --phase PHASE --cell-type TYPE",
    "          --adjust none|batch|age|batch+age --out CSV",
    "Optional: --covariates CSV (required unless adjust=none).",
    sep = "\n")
  if (length(args) == 1L && args[1L] == "--help") {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  if (!length(args) || length(args) %% 2L) stop_clean(usage)
  keys <- args[seq(1L, length(args), by = 2L)]
  values <- args[seq(2L, length(args), by = 2L)]
  valid <- c("--pseudobulk", "--study", "--phase", "--cell-type",
             "--adjust", "--covariates", "--out")
  required_keys <- setdiff(valid, "--covariates")
  if (anyDuplicated(keys) || !all(keys %in% valid) ||
      !all(required_keys %in% keys)) stop_clean(usage)
  opts <- setNames(values, keys)
  if (any(!nzchar(opts[required_keys]))) stop_clean(usage)
  output <- opts[["--out"]]
  if (file.exists(output) || !dir.exists(dirname(output))) {
    stop_clean("Output exists or output directory does not exist; no overwrites")
  }
  s <- read_stratum(opts[["--pseudobulk"]], opts[["--study"]],
                    opts[["--phase"]], opts[["--cell-type"]])
  covariates <- if ("--covariates" %in% names(opts)) opts[["--covariates"]] else NULL
  results <- run_ql(s, covariate_csv = covariates, adjust = opts[["--adjust"]])
  temporary <- tempfile(pattern = ".edger-", tmpdir = dirname(output), fileext = ".csv")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(results, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, output)) stop_clean("Cannot atomically write model results")
  message("EXPLORATORY ONLY: ", sum(results$status == "tested"), " tested genes; ",
          unique(results$n_case), " case and ", unique(results$n_control),
          " control donors. NOT a clinical or biological validation.")
  invisible(output)
}

if (sys.nframe() == 0L) main()
