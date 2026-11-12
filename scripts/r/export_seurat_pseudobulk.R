#!/usr/bin/env Rscript
# Offline-only Seurat/Matrix adapter for an already authorized LOCAL RDS.
# No internet access, publication, or patient-data sharing is performed.
args <- commandArgs(trailingOnly = TRUE)
args <- sub("^--", "", args)
if (!length(args) || !all(grepl("^[A-Za-z0-9-]+=", args))) {
  stop(paste(
    "Usage: Rscript scripts/r/export_seurat_pseudobulk.R",
    "--rds=/approved/atlas.rds --mapping=/approved/donor_map.csv",
    "--out=data/processed/pseudobulk.csv",
    "--audit=data/processed/export_audit.csv",
    "[--assay=RNA] [--sample-col=sample]",
    "[--celltype-col=AnnotationMain] [--min-cells=20]"
  ))
}
kv <- strsplit(args, "=", fixed = TRUE)
keys <- vapply(kv, function(x) x[[1]], character(1))
values <- vapply(kv, function(x) paste(x[-1], collapse = "="), character(1))
if (anyDuplicated(keys) || any(!nzchar(values))) stop("Duplicate or empty arguments")
opt <- stats::setNames(as.list(values), keys)
required <- c("rds", "mapping", "out", "audit")
if (!all(required %in% names(opt))) {
  stop(paste("Missing required arguments:", paste(setdiff(required, names(opt)), collapse = ",")))
}
if (!all(names(opt) %in% c(required, "assay", "sample-col", "celltype-col", "min-cells"))) {
  stop("Unknown argument; please use documented options only")
}
opt_get <- function(name, default) if (is.null(opt[[name]])) default else opt[[name]]
if (!file.exists(opt[["rds"]]) || !file.exists(opt[["mapping"]])) {
  stop("Local RDS and audited donor mapping CSV are required")
}
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("Please install the supported version of SeuratObject in your approved local R environment")
}
cmd <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(cmd)) stop("Run this file with Rscript")
script_path <- sub("^--file=", "", cmd[[1]])
source(file.path(dirname(script_path), "pseudobulk_core.R"), local = TRUE)

obj <- readRDS(opt[["rds"]])
if (!inherits(obj, "Seurat")) stop("Expected an audited Seurat object")
assay <- opt_get("assay", "RNA")
if (!(assay %in% names(obj@assays))) stop("Selected raw assay is absent")
counts <- tryCatch(
  SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
  error = function(err) {
    stop(paste("Could not extract one raw counts layer:",
               conditionMessage(err),
               "Verify Seurat version/layers manually; never silently substitute normalized data."))
  }
)
if (nrow(counts) == 0L || ncol(counts) == 0L) {
  stop("Raw counts layer is empty")
}
meta <- obj[[]]
map <- utils::read.csv(opt[["mapping"]], stringsAsFactors = FALSE, check.names = FALSE)
min_cells <- as.integer(opt_get("min-cells", "20"))
od_export_pseudobulk(
  counts = counts, cell_metadata = meta, sample_map = map,
  output_csv = opt[["out"]], audit_csv = opt[["audit"]],
  sample_col = opt_get("sample-col", "sample"),
  celltype_col = opt_get("celltype-col", "AnnotationMain"),
  min_cells = min_cells
)
