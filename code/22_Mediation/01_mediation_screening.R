#!/bin/env Rscript

## Louise Huuki-Myers & Bernie Mulvey, June 2025
## adapting for mediation screening by Geo Pertea, January 2026
##
## We implement 3 mediation screening scenarios for the carrier status -> Oligo 3. gene expression regression:
## what changes is the gene expression data used in each scenario:
# 1. APOE-associated LC domain DEGs (NM+/- and NM intensity DEGs)
# 2. APOE-associated DEGs from the LC astrocyte spatial domain analysis
# 3. APOE-associated ERC Astrocyte DEGs

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("SummarizedExperiment")
library("BiocParallel")
library("here")
library("getopt")
library("sessioninfo")
library("qs2")


FDRthr=0.05
## FDRthr=0.1 ## loosen for screening

medSelFDR <- 0.05
drop_donors_global <- c("Br1289")
drop_donors_lc <- character(0)

# Helpers
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) < 1 || all(is.na(x))) y else x
}

md5_string <- function(x) {
    tf <- tempfile()
    writeLines(x, con = tf, useBytes = TRUE)
    md5 <- as.character(tools::md5sum(tf))
    unlink(tf)
    md5
}

file_fingerprint <- function(path) {
    info <- file.info(path)
    if (nrow(info) < 1 || is.na(info$size[1]) || is.na(info$mtime[1])) {
        return(sprintf("%s|missing", normalizePath(path, winslash = "/", mustWork = FALSE)))
    }
    sprintf(
        "%s|%s|%s",
        normalizePath(path, winslash = "/", mustWork = FALSE),
        as.character(info$size[1]),
        format(info$mtime[1], tz = "UTC", usetz = TRUE)
    )
}

hash_dt <- function(dt, cols = names(dt)) {
    tf <- tempfile(fileext = ".rds")
    on.exit(unlink(tf), add = TRUE)
    saveRDS(as.data.frame(dt[, ..cols]), file = tf, version = 2)
    as.character(tools::md5sum(tf))
}

hash_matrix <- function(mat) {
    tf <- tempfile(fileext = ".rds")
    on.exit(unlink(tf), add = TRUE)
    saveRDS(
        list(
            rownames = rownames(mat),
            colnames = colnames(mat),
            values = as.vector(mat)
        ),
        file = tf,
        version = 2
    )
    as.character(tools::md5sum(tf))
}

make_donor_key <- function(brnums) {
    donors <- sort(unique(as.character(brnums)))
    donor_str <- paste(donors, collapse = ",")
    donor_key <- md5_string(donor_str)
    list(donors = donors, donor_str = donor_str, donor_key = donor_key)
}

baseline_formula_str <- "0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio"
mediation_formula_str <- paste(baseline_formula_str, "+ med_vec")
baseline_formula <- as.formula(paste("~", baseline_formula_str))
mediation_formula <- as.formula(paste("~", mediation_formula_str))

# Use single baseline (max donors) for all comparisons instead of separate baselines per donor set
forceSingleBaseline <- TRUE

build_design <- function(sample_df, formula_obj, batch_col, med_vec = NULL) {
    sample_df <- as.data.frame(sample_df)
    if ("APOE_syn" %in% names(sample_df)) sample_df$APOE_syn <- droplevels(factor(sample_df$APOE_syn))
    if ("Sex" %in% names(sample_df)) sample_df$Sex <- droplevels(factor(sample_df$Sex))
    if (!(batch_col %in% names(sample_df))) {
        stop(sprintf("Batch column '%s' not found in samples.", batch_col))
    }
    sample_df[[batch_col]] <- droplevels(factor(sample_df[[batch_col]]))
    if (!is.null(med_vec)) sample_df$med_vec <- med_vec
    des <- model.matrix(formula_obj, data = sample_df)
    des <- as.data.frame(des)
    colnames(des) <- gsub(colnames(des), pattern = "_syn", replacement = "_")
    list(design = des, sample_df = sample_df)
}

write_run_info <- function(file, info) {
    if (is.null(names(info)) || any(names(info) == "")) {
        stop("write_run_info requires a named list.")
    }
    lines <- vapply(names(info), function(k) {
        v <- info[[k]]
        if (is.null(v) || length(v) < 1) v <- ""
        sprintf("%s: %s", k, paste(as.character(v), collapse = ","))
    }, character(1))
    writeLines(lines, con = file)
}

read_run_info <- function(file) {
    if (!file.exists(file)) return(NULL)
    lines <- readLines(file, warn = FALSE)
    if (length(lines) < 1) return(NULL)
    kv <- strsplit(lines, ": ", fixed = TRUE)
    keys <- vapply(kv, `[[`, character(1), 1)
    vals <- vapply(kv, function(x) paste(x[-1], collapse = ": "), character(1))
    out <- as.list(vals)
    names(out) <- keys
    out$run_label <- as.character(out$run_label %||% out$run %||% "")
    out$donor_str <- as.character(out$donor_str %||% out$donors %||% "")
    out$donors <- if (nzchar(out$donor_str)) strsplit(out$donor_str, ",", fixed = TRUE)[[1]] else character(0)
    out$n_donors <- as.integer(out$n_donors %||% NA_character_)
    out$n_outcome_genes <- as.integer(out$n_outcome_genes %||% NA_character_)
    out$iter <- as.integer(out$iter %||% NA_character_)
    out
}

cache_matches <- function(cached, expected) {
    if (is.null(cached)) return(FALSE)
    for (nm in names(expected)) {
        cached_val <- as.character(cached[[nm]] %||% "")
        expected_val <- as.character(expected[[nm]] %||% "")
        if (!identical(cached_val, expected_val)) return(FALSE)
    }
    TRUE
}

write_fdr_table <- function(dt, file) {
    ## Never emit empty files; remove stale file when no rows.
    if (nrow(dt) < 1) {
        if (file.exists(file)) unlink(file)
        return(FALSE)
    }
    data.table::fwrite(dt, file = file, sep = "\t", compress = "gzip")
    TRUE
}

write_nonempty_lines <- function(lines, file) {
    lines <- as.character(lines)
    lines <- lines[!is.na(lines) & nzchar(lines)]
    if (length(lines) < 1) {
        if (file.exists(file)) unlink(file)
        return(FALSE)
    }
    writeLines(lines, con = file)
    TRUE
}

# Import command-line parameters
medopt <- matrix(
    c("mediation", "m", "1", "character", "Mediation type: erc_astro, lc_astro, lc_nm"),
    ncol = 5, byrow = TRUE
)
opt <- tryCatch(getopt(medopt), error = function(e) list())
if (is.null(opt$mediation) || is.na(opt$mediation)) opt$mediation <- "lc_astro"

if (!(opt$mediation %in% c("erc_astro", "lc_astro", "lc_nm"))) {
    stop("Invalid mediation type. Choose from: erc_astro, lc_astro, lc_nm")
}
if (opt$mediation == "lc_nm") {
    stop(sprintf("Mediation type '%s' is not implemented yet.", opt$mediation))
}

#workers <- as.integer(Sys.getenv("MEDIATION_WORKERS", "10"))
#if (is.na(workers) || workers < 1) workers <- 1
workers <- as.integer(Sys.getenv("MEDIATION_WORKERS", "10"))
if (is.na(workers) || workers < 1) workers <- 10L

parse_int_csv <- function(x) {
    x <- trimws(as.character(x %||% ""))
    if (!nzchar(x)) return(integer(0))
    toks <- strsplit(x, ",", fixed = TRUE)[[1]]
    toks <- trimws(toks)
    toks <- toks[nzchar(toks)]
    if (length(toks) < 1) return(integer(0))
    vals <- suppressWarnings(as.integer(toks))
    if (any(is.na(vals))) {
        stop(sprintf("Invalid MEDIATION_TEST_ROWS value: %s", x))
    }
    unique(vals)
}

test_rows <- parse_int_csv(Sys.getenv("MEDIATION_TEST_ROWS", ""))
#test_rows <- c(1,10,20) ## alternative manual test rows
if (length(test_rows) > 0) {
        message("Mediation test rows provided: ", paste(test_rows, collapse = ", "))
}

message(Sys.time(), sprintf(" - Mediation type = %s (workers = %i)", opt$mediation, workers))

#### Paths ####
ddir <- here("processed-data")
mdir <- here("processed-data", "22_Mediation")

pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn", "sce_pseudo_DGE-cell_type_anno.RDS")

out_dir_env <- Sys.getenv("MEDIATION_OUTDIR", "")
out_dir <- if (nzchar(trimws(out_dir_env))) {
    out_dir_env
} else {
    here("code", "22_Mediation", sprintf("out-%s", opt$mediation))
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
baseline_dir <- file.path(out_dir, "baseline")
dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
run_input_files <- c(pb_fn)

#### Load X->Y pseudo bulk data ####
sce_pb <- readRDS(pb_fn)
## always drop Br1289
sce_pb <- sce_pb[, !as.character(sce_pb$BrNum) %in% drop_donors_global]

#### Outcome cluster setup (same for every scenario) ####
## we can use dge_base colData() for all donor metadata

out_clus <- "Oligo.3"
batch <- "exp_round"
sce_base <- sce_pb[, sce_pb$registration_variable == out_clus]
donor_meta <- as.data.frame(colData(sce_base))


med_degs <- NULL ## placeholder for the selected mediator DEGs, gene metadata
med_expr <- NULL ## placeholder for the mediator gene expression data (logcounts)
## these should be loaded depending on mediation type

if (opt$mediation == "erc_astro") { ## "green" design 3: Mediator = ERC Astro DEGs
  degs_fn <- file.path(ddir, "13_compile_DGE/01_compile_DGE/sn_fine", "DGE_results_carrier_sn_fine.Rds")
  run_input_files <- c(run_input_files, degs_fn)
  pdge <- readRDS(degs_fn)
  med_degs <- as.data.table(pdge)[
      grepl("^Astro", cluster) & vlmf_adj.P.Val < medSelFDR & contrast == "carrier",
      .(cluster, gene_id, gene_name, gene_type, vlmf_logFC, vlmf_AveExpr, vlmf_t, vlmf_P.Value, vlmf_adj.P.Val, cell_type_broad)
  ]
  setorder(med_degs, cluster, gene_id)
  med_degs[, med_gene_id := paste0(cluster, "|", gene_id)] ## unique mediator ids, matches med_expr rownames

  ## only for this mediation scenario, the mediator expression is in the same sce_pb object
  ## prepare mediator expression matrix (med_gene_id x BrNum) from sce_pb logcounts
  med_brnums_all <- unique(as.character(colData(sce_pb)$BrNum))
  ## > head(erc_deg_summary)
  #    V1     cell_type_broad  cluster   model   n_genes   nDEGs_FDR05   nUP   nDown
  #   <int>       <char>       <char>   <char>     <int>      <int>     <int>  <int>
  #1:  1         Astro        Astro.1   carrier   16867         7         5     2
  #2:  2         Astro        Astro.2   carrier    9994        58        10    48
  #3:  3         Astro        Astro.3   carrier   11282        23         5    18

  ## Note: for Astro.1 cluster, only 29 out of 30 donors had enough expression data
  med_expr <- matrix(NA_real_, nrow = nrow(med_degs), ncol = length(med_brnums_all))
  rownames(med_expr) <- med_degs$med_gene_id
  colnames(med_expr) <- med_brnums_all

  for (clus in unique(med_degs$cluster)) {
      row_idx <- which(med_degs$cluster == clus)
      gene_ids <- as.character(med_degs$gene_id[row_idx])
      missing_ids <- setdiff(gene_ids, rownames(sce_pb))
      if (length(missing_ids) > 0) {
          stop(sprintf("Missing gene_id(s) in sce_pb for cluster %s: %s", clus, paste(missing_ids, collapse = ", ")))
      }
      clus_sce <- sce_pb[gene_ids, as.character(sce_pb$registration_variable) == clus]
      clus_brnum <- as.character(colData(clus_sce)$BrNum)
      if (anyDuplicated(clus_brnum)) {
          stop(sprintf("Duplicate BrNum entries found for cluster %s in sce_pb.", clus))
      }
      clus_expr <- as.matrix(assay(clus_sce, "logcounts"))
      rownames(clus_expr) <- paste0(clus, "|", rownames(clus_expr))
      med_expr[rownames(clus_expr), clus_brnum] <- clus_expr
  }

  if (length(test_rows) > 0) {
      if (any(test_rows < 1 | test_rows > nrow(med_degs))) {
          stop(sprintf("MEDIATION_test_ROWS has out-of-range indices. Valid range: 1..%d", nrow(med_degs)))
      }
      med_degs <- med_degs[test_rows]
      med_expr <- med_expr[med_degs$med_gene_id, , drop = FALSE]
      message(Sys.time(), sprintf(" - test test enabled: using med_degs rows %s", paste(test_rows, collapse = ",")))
  }
} else if (opt$mediation == "lc_astro") { ## "blue" design 2: Mediator = LC Astro DEGs from LC spatial domains
    ## Load LC Astro DEGs (FDR < 0.05)
    f_lca_dge <- file.path(mdir, "LC_Astro_E4vE2_DGEout_Bernie.rds")
    if (!file.exists(f_lca_dge)) {
        stop(sprintf("LC Astro DEG file not found: %s", f_lca_dge))
    }
    run_input_files <- c(run_input_files, f_lca_dge)
    lca_dge <- as.data.table(readRDS(f_lca_dge))
    med_degs <- lca_dge[adj.P.Val < medSelFDR]
    if (nrow(med_degs) < 1) {
        stop(sprintf("No LC Astro DEGs with FDR < %.4f found.", medSelFDR))
    }
    ## Standardize column names and add required fields
    med_degs[, cluster := "LCAstro"]
    med_degs[, med_gene_id := paste0("LCAstro|", gene_id)]
    setorder(med_degs, gene_id)

    ## Load LC Astro pseudobulk counts (DGEList)
    f_lca_pbcounts <- file.path(mdir, "lc_astro_pseudobulk_counts_by_donor.rds")
    if (!file.exists(f_lca_pbcounts)) {
        stop(sprintf("LC Astro pseudobulk counts file not found: %s", f_lca_pbcounts))
    }
    run_input_files <- c(run_input_files, f_lca_pbcounts)
    dgl_lca <- readRDS(f_lca_pbcounts)
    ## Determine overlapping donors between LC Astro and ERC Oligo.3 outcome
    ## after leave-one-out testing, let's drop BR5529 from dgl_lca
    #drop_donors_lc <- c("Br5529")
    drop_donors_lc <- c("Br5529", "Br6423") ## drop Br6423 too, which has very low coverage in LC Astro and is an outlier in PCA
    dgl_lca <- dgl_lca[, !colnames(dgl_lca) %in% drop_donors_lc]
    #dgl_lca <- dgl_lca[, !colnames(dgl_lca) %in% c("Br5529", "Br6423")]
    lca_brnums <- colnames(dgl_lca)
    erc_brnums <- as.character(donor_meta$BrNum)
    common_donors <- intersect(lca_brnums, erc_brnums)
    if (length(common_donors) < 10) {
        stop(sprintf("Too few overlapping donors between LC Astro and ERC Oligo.3: %d", length(common_donors)))
    }
    message(Sys.time(), sprintf(" - LC-ERC overlapping donors: %d", length(common_donors)))

    ## Subset LC counts to overlapping donors and filter by expression
    dgl_lca <- dgl_lca[, common_donors]
    ## Build design matrix for filtering using subsetted donor_meta
    lca_donor_meta <- donor_meta[donor_meta$BrNum %in% common_donors, ]
    lca_donor_meta <- lca_donor_meta[match(common_donors, lca_donor_meta$BrNum), ]
    lca_design <- model.matrix(baseline_formula, data = lca_donor_meta)

    ## Filter and normalize LC counts
    dgl_lca <- edgeR::calcNormFactors(dgl_lca)
    keep_lca <- edgeR::filterByExpr(dgl_lca, design = lca_design)
    dgl_lca <- dgl_lca[keep_lca, , keep.lib.sizes = FALSE]
    dgl_lca <- edgeR::calcNormFactors(dgl_lca)

    ## compute log-CPM - will be used for the mediation vector
    lca_logcpm <- edgeR::cpm(dgl_lca, log = TRUE)

    ## Check for missing DEGs in filtered counts
    gene_ids <- as.character(med_degs$gene_id)
    missing_ids <- setdiff(gene_ids, rownames(lca_logcpm))
    if (length(missing_ids) > 0) {
        warning(sprintf("Dropping %d LC Astro DEGs not in filtered counts: %s",
                        length(missing_ids), paste(head(missing_ids, 5), collapse = ", ")))
        med_degs <- med_degs[!(gene_id %in% missing_ids)]
        gene_ids <- as.character(med_degs$gene_id)
    }
    if (nrow(med_degs) < 1) {
        stop("No LC Astro DEGs remaining after filtering.")
    }

    ## Build med_expr matrix (med_gene_id x BrNum)
    med_expr <- lca_logcpm[gene_ids, , drop = FALSE]
    rownames(med_expr) <- med_degs$med_gene_id
    ## Apply test_rows if specified
    if (length(test_rows) > 0) {
        if (any(test_rows < 1 | test_rows > nrow(med_degs))) {
            stop(sprintf("test_rows has out-of-range indices. Valid range: 1..%d", nrow(med_degs)))
        }
        med_degs <- med_degs[test_rows]
        med_expr <- med_expr[med_degs$med_gene_id, , drop = FALSE]
        message(Sys.time(),
             sprintf(" - test enabled: using med_degs rows %s", paste(test_rows, collapse = ",")))
    }
} else if (opt$mediation == "lc_nm") { ## "orange" design, with NM data
} else {
  stop(sprintf("Mediation type '%s' not implemented yet.", opt$mediation))
}

if (is.null(med_degs) || is.null(med_expr)) {
    stop(sprintf("med_degs/med_expr not initialized for mediation type '%s'.", opt$mediation))
}
if (nrow(med_degs) < 1) {
    stop(sprintf("No mediator genes found for mediation type '%s'.", opt$mediation))
}
if (!all(c("cluster", "gene_id", "gene_name", "med_gene_id") %in% names(med_degs))) {
    stop("med_degs is missing required columns: cluster, gene_id, gene_name, med_gene_id.")
}

## Canonical mediator ordering and stable iteration ids for file naming.
setorder(med_degs, cluster, gene_id, gene_name)
med_degs[, med_gene_id := paste0(cluster, "|", gene_id)]
if (anyDuplicated(med_degs$med_gene_id)) {
    dup <- unique(med_degs$med_gene_id[duplicated(med_degs$med_gene_id)])
    stop(sprintf("Duplicate mediator IDs found (med_gene_id), first examples: %s", paste(head(dup, 5), collapse = ", ")))
}
med_degs[, run := sprintf("%s|%s|%s", cluster, gene_id, gene_name)]
med_degs[, iter := seq_len(.N)]
if (!all(med_degs$med_gene_id %in% rownames(med_expr))) {
    miss <- setdiff(med_degs$med_gene_id, rownames(med_expr))
    stop(sprintf("Missing mediators in med_expr, first examples: %s", paste(head(miss, 5), collapse = ", ")))
}
med_expr <- med_expr[med_degs$med_gene_id, , drop = FALSE]
med_expr <- med_expr[, sort(colnames(med_expr)), drop = FALSE]
message(Sys.time(), sprintf(" - Mediators selected for %s: %i", opt$mediation, nrow(med_degs)))
message(Sys.time(), sprintf(" - Mediator donor columns: %i", ncol(med_expr)))

## subset dge_base to the column names in med_expr -- the actual donors to keep
dge_base <- sce_base[, sce_base$BrNum %in% colnames(med_expr)]

## Filter low expression genes once (baseline design)
fn_dge_base <- file.path(out_dir, "dge_base.qs2")
dge_base <- edgeR::calcNormFactors(dge_base)
dge_base_filter_design <- model.matrix(baseline_formula, data = as.data.frame(dge_base$samples))
keep <- edgeR::filterByExpr(dge_base, design = dge_base_filter_design)
dge_base <- dge_base[keep, , keep.lib.sizes = FALSE]
dge_base <- edgeR::calcNormFactors(dge_base)
outcome_donor_info <- make_donor_key(as.character(dge_base$samples$BrNum))
outcome_gene_hash <- md5_string(paste(rownames(dge_base), collapse = "\n"))
med_degs_hash <- hash_dt(med_degs, c("iter", "cluster", "gene_id", "gene_name", "med_gene_id", "run"))
med_expr_hash <- hash_matrix(med_expr)
run_signature_lines <- c(
    sprintf("mediation=%s", opt$mediation),
    sprintf("FDRthr=%s", formatC(FDRthr, digits = 8, format = "fg")),
    sprintf("medSelFDR=%s", formatC(medSelFDR, digits = 8, format = "fg")),
    sprintf("forceSingleBaseline=%s", forceSingleBaseline),
    sprintf("batch=%s", batch),
    sprintf("baseline_formula=%s", baseline_formula_str),
    sprintf("mediation_formula=%s", mediation_formula_str),
    sprintf("drop_donors_global=%s", paste(drop_donors_global, collapse = ",")),
    sprintf("drop_donors_lc=%s", paste(drop_donors_lc, collapse = ",")),
    sprintf("test_rows=%s", paste(test_rows, collapse = ",")),
    sprintf("n_mediators=%d", nrow(med_degs)),
    sprintf("med_degs_hash=%s", med_degs_hash),
    sprintf("med_expr_hash=%s", med_expr_hash),
    sprintf("outcome_gene_hash=%s", outcome_gene_hash),
    sprintf("outcome_donor_key=%s", outcome_donor_info$donor_key),
    sprintf("outcome_donor_str=%s", outcome_donor_info$donor_str),
    paste(sprintf("input_file=%s", vapply(run_input_files, file_fingerprint, character(1))), collapse = "\n")
)
run_hash <- md5_string(paste(run_signature_lines, collapse = "\n"))
run_started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
message(Sys.time(), sprintf(" - Run hash: %s", run_hash))

## save run-state objects; always overwrite latest to avoid stale downstream reads.
qs_save(med_degs, file.path(out_dir, "med_degs.qs2"))
qs_save(med_expr, file.path(out_dir, "med_expr.qs2"))
qs_save(dge_base, fn_dge_base)
qs_save(med_degs, file.path(out_dir, sprintf("med_degs_%s.qs2", run_hash)))
qs_save(med_expr, file.path(out_dir, sprintf("med_expr_%s.qs2", run_hash)))
qs_save(dge_base, file.path(out_dir, sprintf("dge_base_%s.qs2", run_hash)))

writeLines(run_signature_lines, con = file.path(out_dir, "run_signature_latest.txt"))
writeLines(run_signature_lines, con = file.path(out_dir, sprintf("run_signature_%s.txt", run_hash)))
data.table::fwrite(
    med_degs[, .(iter, run, cluster, gene_id, gene_name, med_gene_id)],
    file = file.path(out_dir, "mediation_mediator_index.tsv"),
    sep = "\t"
)
data.table::fwrite(
    med_degs[, .(iter, run, cluster, gene_id, gene_name, med_gene_id)],
    file = file.path(out_dir, sprintf("mediation_mediator_index.%s.tsv", run_hash)),
    sep = "\t"
)

## outcome gene mapping for per-run list outputs
get_out_gene_map <- function(dge_obj) {
    gene_name_vec <- NULL
    if (inherits(dge_obj, "DGEList")) {
        if (!is.null(dge_obj$genes) && "gene_name" %in% colnames(dge_obj$genes)) {
            gene_name_vec <- as.character(dge_obj$genes$gene_name)
        }
    } else if ("gene_name" %in% colnames(rowData(dge_obj))) {
        gene_name_vec <- as.character(rowData(dge_obj)$gene_name)
    }
    if (is.null(gene_name_vec)) {
        gene_name_vec <- rep(NA_character_, nrow(dge_obj))
    }
    out_gene_map <- data.table::data.table(
        outcome_gene_id = rownames(dge_obj),
        outcome_gene_name = gene_name_vec
    )
    out_gene_map <- unique(out_gene_map, by = "outcome_gene_id")
    data.table::setkey(out_gene_map, outcome_gene_id)
    out_gene_map
}
out_gene_map <- get_out_gene_map(dge_base)

## precompute per-mediator donor plans from med_expr (before running models)
out_brnum_all <- as.character(dge_base$samples$BrNum)
mediator_plan_dt <- med_degs[, {
    if (!(med_gene_id %in% rownames(med_expr))) {
        stop(sprintf("med_gene_id missing in med_expr during plan build: %s", med_gene_id))
    }
    med_vals <- as.numeric(med_expr[med_gene_id, out_brnum_all])
    donors_i <- out_brnum_all[!is.na(med_vals)]
    donor_info <- make_donor_key(donors_i)
    list(
        donors = list(donor_info$donors),
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str,
        n_donors = length(donor_info$donors),
        eligible = length(donor_info$donors) >= 10L
    )
}, by = .(iter, run_label = run, med_gene_id, gene_id, gene_name)]
setorder(mediator_plan_dt, iter)

n_skipped_plan <- sum(!mediator_plan_dt$eligible)
if (n_skipped_plan > 0) {
    message(Sys.time(), sprintf(" - Skipping %d mediators with <10 donors (pre-check)", n_skipped_plan))
}
mediator_plan_run <- mediator_plan_dt[eligible == TRUE]
if (nrow(mediator_plan_run) < 1) {
    stop("No eligible mediators remain after donor-overlap pre-check.")
}

run_baseline <- function(donor_key, donors, donor_str) {
    out_file <- file.path(baseline_dir, sprintf("baseline_%s.tab.gz", donor_key))
    sufdr <- stringr::str_split_1(as.character(FDRthr), "\\.")[[2]]
    out_fdr_file <- file.path(baseline_dir, sprintf("baseline_%s.fdr%s.tab.gz", donor_key, sufdr))
    baseline_info_file <- file.path(baseline_dir, sprintf("baseline_%s.txt", donor_key))

    out_brnum <- as.character(dge_base$samples$BrNum)
    common_donors <- intersect(out_brnum, donors)
    if (length(common_donors) < 10) {
        warning(sprintf("Too few common donors (%d) for baseline %s", length(common_donors), donor_key))
        return(NULL)
    }
    donor_info <- make_donor_key(common_donors)
    baseline_hash <- md5_string(paste(
        c(
            run_hash,
            "baseline",
            donor_info$donor_key,
            donor_info$donor_str,
            baseline_formula_str,
            batch,
            outcome_gene_hash,
            as.character(nrow(dge_base)),
            as.character(ncol(dge_base))
        ),
        collapse = "\n"
    ))
    expected_cache <- list(
        run_hash = run_hash,
        baseline_hash = baseline_hash,
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str
    )
    if (file.exists(out_file) && file.exists(baseline_info_file)) {
        cached <- read_run_info(baseline_info_file)
        if (cache_matches(cached, expected_cache)) {
            message(Sys.time(), sprintf(" - Using cached baseline: %s", out_file))
            tt <- data.table::fread(out_file)
            baseline_sig <- tt[adj.P.Val < FDRthr, outcome_gene_id]
            if (!file.exists(out_fdr_file)) {
                write_fdr_table(tt[adj.P.Val < FDRthr], out_fdr_file)
            }
            return(list(
                donor_key = donor_info$donor_key,
                donor_str = donor_info$donor_str,
                n_donors = length(common_donors),
                baseline_hash = baseline_hash,
                out_file = out_file,
                sig_genes = list(baseline_sig)
            ))
        }
        message(Sys.time(), sprintf(" - Baseline cache invalidated, recomputing: %s", out_file))
    }

    dge_sub <- dge_base[, out_brnum %in% common_donors]
    design_obj <- build_design(dge_sub$samples, baseline_formula, batch)
    des <- design_obj$design
    sample_df <- design_obj$sample_df

    v.swt <- voomLmFit(
        dge_sub,
        design = des,
        block = as.factor(sample_df[[batch]]),
        adaptive.span = TRUE,
        sample.weights = TRUE
    )

    cont <- makeContrasts(
        carrier = "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
        levels = des
    )

    v.swt.fit <- contrasts.fit(v.swt, contrasts = cont)
    v.swt.fit.e <- eBayes(v.swt.fit)

    tt <- topTable(v.swt.fit.e, coef = "carrier", number = Inf, adjust.method = "BH")
    tt <- as.data.table(tt, keep.rownames = "outcome_gene_id")
    tt[, `:=`(run = sprintf("baseline|%s", donor_key), cluster = out_clus, contrast = "carrier",
              n_donors = length(common_donors), donor_key = donor_key)]
    setcolorder(tt, c("run", "cluster", "contrast", "outcome_gene_id", "n_donors", "donor_key",
                      setdiff(names(tt), c("run", "cluster", "contrast", "outcome_gene_id", "n_donors", "donor_key"))))

    data.table::fwrite(tt, file = out_file, sep = "\t", compress = "gzip")
    write_fdr_table(tt[adj.P.Val < FDRthr], out_fdr_file)
    write_run_info(baseline_info_file, list(
        run_hash = run_hash,
        baseline_hash = baseline_hash,
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str,
        donors = donor_info$donor_str,
        n_donors = length(common_donors),
        n_outcome_genes = nrow(dge_sub),
        design_formula = baseline_formula_str,
        design_dim = sprintf("%dx%d", nrow(des), ncol(des))
    ))

    baseline_sig <- tt[adj.P.Val < FDRthr, outcome_gene_id]
    list(
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str,
        n_donors = length(common_donors),
        baseline_hash = baseline_hash,
        out_file = out_file,
        sig_genes = list(baseline_sig)
    )
}

## compute baseline(s) before run_one parallel phase
bpparam <- BiocParallel::MulticoreParam(workers = workers)
donor_sets <- unique(mediator_plan_run[, .(donor_key, donor_str, donors)], by = "donor_key")
donor_sets[, n_donors := lengths(donors)]
setorder(donor_sets, -n_donors, donor_key)

if (forceSingleBaseline) {
    donor_sets_to_run <- donor_sets[1]
    single_baseline_key <- donor_sets_to_run$donor_key
    message(Sys.time(), sprintf(" - Using single baseline with %d donors (forceSingleBaseline=TRUE)",
                                 donor_sets_to_run$n_donors[1]))
} else {
    donor_sets_to_run <- donor_sets
    single_baseline_key <- NULL
}

baseline_results <- BiocParallel::bplapply(
    seq_len(nrow(donor_sets_to_run)),
    function(i) run_baseline(donor_sets_to_run$donor_key[i], donor_sets_to_run$donors[[i]], donor_sets_to_run$donor_str[i]),
    BPPARAM = bpparam
)
baseline_results <- Filter(Negate(is.null), baseline_results)
if (length(baseline_results) < 1) {
    stop("No baseline runs completed successfully.")
}

baseline_dt <- data.table::rbindlist(baseline_results, fill = TRUE)
baseline_map <- setNames(lapply(baseline_results, function(x) x$sig_genes[[1]]),
                         vapply(baseline_results, `[[`, character(1), "donor_key"))

if (forceSingleBaseline && !is.null(single_baseline_key)) {
    single_baseline_genes <- baseline_map[[single_baseline_key]]
    for (dk in mediator_plan_run$donor_key) {
        baseline_map[[dk]] <- single_baseline_genes
    }
}
missing_baselines <- setdiff(unique(mediator_plan_run$donor_key), names(baseline_map))
if (length(missing_baselines) > 0) {
    stop(sprintf("Missing baseline mappings for donor keys: %s", paste(head(missing_baselines, 5), collapse = ", ")))
}

baseline_index_dt <- baseline_dt[, .(run_hash = run_hash, donor_key, donor_str, n_donors, baseline_hash, out_file)]
data.table::fwrite(
    baseline_index_dt,
    file = file.path(out_dir, "baseline_index.tsv"),
    sep = "\t"
)
data.table::fwrite(
    baseline_index_dt,
    file = file.path(out_dir, sprintf("baseline_index.%s.tsv", run_hash)),
    sep = "\t"
)

summarize_vs_baseline <- function(tt, ttM, run_label, donor_key) {
    carrier_sig <- tt[adj.P.Val < FDRthr, outcome_gene_id]
    med_sig <- ttM[adj.P.Val < FDRthr, outcome_gene_id]
    base_genes <- baseline_map[[donor_key]]
    if (is.null(base_genes)) {
        stop(sprintf("Missing baseline for donor_key %s", donor_key))
    }
    lost <- setdiff(base_genes, carrier_sig)
    gained <- setdiff(carrier_sig, base_genes)
    attenuated <- lost
    step3a_genes <- intersect(med_sig, base_genes)
    mediated <- intersect(attenuated, med_sig)

    lost_dt <- out_gene_map[J(lost), on = "outcome_gene_id"]
    gained_dt <- out_gene_map[J(gained), on = "outcome_gene_id"]
    mediated_dt <- out_gene_map[J(mediated), on = "outcome_gene_id"]
    lost_labels <- ifelse(is.na(lost_dt$outcome_gene_name) | lost_dt$outcome_gene_name == "",
                          lost_dt$outcome_gene_id,
                          paste0(lost_dt$outcome_gene_id, "|", lost_dt$outcome_gene_name))
    gained_labels <- ifelse(is.na(gained_dt$outcome_gene_name) | gained_dt$outcome_gene_name == "",
                            gained_dt$outcome_gene_id,
                            paste0(gained_dt$outcome_gene_id, "|", gained_dt$outcome_gene_name))
    mediated_labels <- ifelse(is.na(mediated_dt$outcome_gene_name) | mediated_dt$outcome_gene_name == "",
                              mediated_dt$outcome_gene_id,
                              paste0(mediated_dt$outcome_gene_id, "|", mediated_dt$outcome_gene_name))

    file_label <- gsub("\\|", "_", run_label)
    loss_file <- file.path(out_dir, sprintf("%s.loss.txt", file_label))
    gain_file <- file.path(out_dir, sprintf("%s.gain.txt", file_label))
    mediated_file <- file.path(out_dir, sprintf("%s.mediated.lst", file_label))

    wrote_loss <- write_nonempty_lines(sort(as.character(lost_labels)), loss_file)
    wrote_gain <- write_nonempty_lines(sort(as.character(gained_labels)), gain_file)
    wrote_mediated <- write_nonempty_lines(sort(as.character(mediated_labels)), mediated_file)

    list(
        baselineDEGs = length(base_genes),
        mediationCarrierDEGs = length(carrier_sig),
        mediationMDEGs = length(med_sig),
        medSigBaselineDEGs = length(step3a_genes),
        mediated_n = length(mediated),
        pass_both = length(mediated) > 0,
        loss = length(lost),
        gain = length(gained),
        loss_file = if (wrote_loss) loss_file else NA_character_,
        gain_file = if (wrote_gain) gain_file else NA_character_,
        mediated_file = if (wrote_mediated) mediated_file else NA_character_
    )
}

run_one <- function(i) {
    plan_row <- mediator_plan_run[i]
    iter_id <- as.integer(plan_row$iter[[1]])

    ## Use a non-column symbol on RHS to avoid data.table scoping collisions.
    med_row <- med_degs[iter == iter_id]
    if (nrow(med_row) != 1) {
        stop(sprintf("Could not resolve exactly one mediator row for iter %d (found %d).", iter_id, nrow(med_row)))
    }

    run_label <- as.character(med_row$run)
    med_gene_id <- as.character(med_row$med_gene_id)
    gene_id <- as.character(med_row$gene_id)
    gene_name <- as.character(med_row$gene_name)

    expected_run_label <- sprintf("%s|%s|%s", med_row$cluster, med_row$gene_id, med_row$gene_name)
    expected_med_gene_id <- sprintf("%s|%s", med_row$cluster, med_row$gene_id)
    if (!identical(run_label, expected_run_label) || !identical(run_label, as.character(plan_row$run_label))) {
        stop(sprintf("run_label mismatch at iter %d", iter_id))
    }
    if (!identical(med_gene_id, expected_med_gene_id) || !identical(med_gene_id, as.character(plan_row$med_gene_id))) {
        stop(sprintf("med_gene_id mismatch at iter %d", iter_id))
    }
    if (!identical(gene_id, as.character(plan_row$gene_id)) || !identical(gene_name, as.character(plan_row$gene_name))) {
        stop(sprintf("Mediator gene metadata mismatch at iter %d", iter_id))
    }

    donors_i <- as.character(plan_row$donors[[1]])
    if (length(donors_i) < 10L) {
        warning(sprintf("Too few planned donors (%d) for %s", length(donors_i), run_label))
        return(NULL)
    }
    if (!all(donors_i %in% out_brnum_all)) {
        miss <- setdiff(donors_i, out_brnum_all)
        stop(sprintf("Planned donors missing in dge_base for %s: %s", run_label, paste(head(miss, 5), collapse = ", ")))
    }

    donor_info <- make_donor_key(donors_i)
    if (!identical(donor_info$donor_key, as.character(plan_row$donor_key)) ||
        !identical(donor_info$donor_str, as.character(plan_row$donor_str))) {
        stop(sprintf("Donor key/string mismatch at iter %d", iter_id))
    }

    iter <- iter_id
    file_label <- gsub("\\|", "_", run_label)
    out_file <- file.path(out_dir, sprintf("%03d_%s.tab.gz", iter, file_label))
    out_Mfile <- file.path(out_dir, sprintf("%03d_%s.M.tab.gz", iter, file_label))
    sufdr <- stringr::str_split_1(as.character(FDRthr), "\\.")[[2]]
    out_fdr_file <- sub("\\.tab\\.gz$", sprintf(".fdr%s.tab.gz", sufdr), out_file)
    out_Mfdr_file <- sub("\\.tab\\.gz$", sprintf(".fdr%s.tab.gz", sufdr), out_Mfile)
    run_info_file <- file.path(out_dir, sprintf("%03d_%s.txt", iter, file_label))

    dge_sub <- dge_base[, as.character(dge_base$samples$BrNum) %in% donors_i]
    dge_sub <- dge_sub[, match(donors_i, as.character(dge_sub$samples$BrNum))]
    if (!identical(as.character(dge_sub$samples$BrNum), donors_i)) {
        stop(sprintf("Donor ordering mismatch in dge_sub for %s", run_label))
    }

    med_vec <- as.numeric(med_expr[med_gene_id, donors_i])
    if (anyNA(med_vec)) {
        stop(sprintf("Unexpected NA mediator values after donor planning for %s", run_label))
    }

    med_vec_hash <- md5_string(paste(sprintf("%.10f", med_vec), collapse = ","))
    analysis_hash <- md5_string(paste(
        c(
            run_hash,
            run_label,
            med_gene_id,
            donor_info$donor_key,
            donor_info$donor_str,
            as.character(nrow(dge_sub)),
            as.character(ncol(dge_sub)),
            outcome_gene_hash,
            med_vec_hash
        ),
        collapse = "\n"
    ))

    if (file.exists(out_file) && file.exists(out_Mfile) && file.exists(run_info_file)) {
        cached <- read_run_info(run_info_file)
        expected_cache <- list(
            run_hash = run_hash,
            analysis_hash = analysis_hash,
            run_label = run_label,
            med_gene_id = med_gene_id,
            donor_key = donor_info$donor_key,
            donor_str = donor_info$donor_str,
            iter = iter,
            n_outcome_genes = nrow(dge_sub)
        )
        if (cache_matches(cached, expected_cache)) {
            tt <- data.table::fread(out_file)
            ttM <- data.table::fread(out_Mfile)
            if (!file.exists(out_fdr_file)) write_fdr_table(tt[adj.P.Val < FDRthr], out_fdr_file)
            if (!file.exists(out_Mfdr_file)) write_fdr_table(ttM[adj.P.Val < FDRthr], out_Mfdr_file)
            summary_fields <- summarize_vs_baseline(tt, ttM, run_label, donor_info$donor_key)
            message(Sys.time(), sprintf(" - Using cached run: %s", out_file))
            return(c(
                list(
                    iter = iter,
                    out_file = out_file,
                    out_Mfile = out_Mfile,
                    run_label = run_label,
                    med_gene_id = med_gene_id,
                    gene_id = gene_id,
                    gene_name = gene_name,
                    donor_key = donor_info$donor_key,
                    donor_str = donor_info$donor_str,
                    donors = list(donor_info$donors),
                    n_donors = length(donors_i),
                    n_outcome_genes = nrow(dge_sub),
                    run_hash = run_hash,
                    analysis_hash = analysis_hash
                ),
                summary_fields
            ))
        }
        message(Sys.time(), sprintf(" - Cache invalidated, recomputing: %s", out_file))
    }

    design_obj <- build_design(dge_sub$samples, mediation_formula, batch, med_vec = med_vec)
    des <- design_obj$design
    sample_df <- design_obj$sample_df

    v.swt <- voomLmFit(
        dge_sub,
        design = des,
        block = as.factor(sample_df[[batch]]),
        adaptive.span = TRUE,
        sample.weights = TRUE
    )

    cont <- makeContrasts(
        carrier = "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
        med_vec = med_vec,
        levels = des
    )

    v.swt.fit <- contrasts.fit(v.swt, contrasts = cont)
    v.swt.fit.e <- eBayes(v.swt.fit)

    tt <- topTable(v.swt.fit.e, coef = "carrier", number = Inf, adjust.method = "BH")
    tt <- as.data.table(tt, keep.rownames = "outcome_gene_id")
    tt[, `:=`(run = run_label, cluster = out_clus, contrast = "carrier",
              n_donors = length(donors_i))]
    setcolorder(tt, c("run", "cluster", "contrast", "outcome_gene_id", "n_donors",
                      setdiff(names(tt), c("run", "cluster", "contrast", "outcome_gene_id", "n_donors"))))

    ttM <- topTable(v.swt.fit.e, coef = "med_vec", number = Inf, adjust.method = "BH")
    ttM <- as.data.table(ttM, keep.rownames = "outcome_gene_id")
    ttM[, `:=`(run = run_label, cluster = out_clus, contrast = "med_vec",
               n_donors = length(donors_i))]
    setcolorder(ttM, c("run", "cluster", "contrast", "outcome_gene_id", "n_donors",
                       setdiff(names(ttM), c("run", "cluster", "contrast", "outcome_gene_id", "n_donors"))))

    data.table::fwrite(tt, file = out_file, sep = "\t", compress = "gzip")
    write_fdr_table(tt[adj.P.Val < FDRthr], out_fdr_file)
    data.table::fwrite(ttM, file = out_Mfile, sep = "\t", compress = "gzip")
    write_fdr_table(ttM[adj.P.Val < FDRthr], out_Mfdr_file)

    write_run_info(run_info_file, list(
        run_hash = run_hash,
        analysis_hash = analysis_hash,
        run = run_label,
        run_label = run_label,
        iter = iter,
        med_gene_id = med_gene_id,
        gene_id = gene_id,
        gene_name = gene_name,
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str,
        donors = donor_info$donor_str,
        n_donors = length(donors_i),
        n_outcome_genes = nrow(dge_sub),
        design_formula = mediation_formula_str,
        design_dim = sprintf("%dx%d", nrow(des), ncol(des))
    ))

    summary_fields <- summarize_vs_baseline(tt, ttM, run_label, donor_info$donor_key)
    c(
        list(
            iter = iter,
            out_file = out_file,
            out_Mfile = out_Mfile,
            run_label = run_label,
            med_gene_id = med_gene_id,
            gene_id = gene_id,
            gene_name = gene_name,
            donor_key = donor_info$donor_key,
            donor_str = donor_info$donor_str,
            donors = list(donor_info$donors),
            n_donors = length(donors_i),
            n_outcome_genes = nrow(dge_sub),
            run_hash = run_hash,
            analysis_hash = analysis_hash
        ),
        summary_fields
    )
}

message(Sys.time(), sprintf(" - Starting %d workers for %d %s mediation runs",
                            workers, nrow(mediator_plan_run), opt$mediation))
run_results <- BiocParallel::bplapply(seq_len(nrow(mediator_plan_run)), run_one, BPPARAM = bpparam)
run_results <- Filter(Negate(is.null), run_results)
message(Sys.time(), sprintf(" - Runs completed (%i files)", length(run_results)))

run_info_dt <- data.table::rbindlist(run_results, fill = TRUE)
setorder(run_info_dt, iter)
if (nrow(run_info_dt) < 1) {
    stop("No mediation runs completed successfully.")
}

data.table::fwrite(
    run_info_dt[, .(
        iter, run_label, med_gene_id, gene_id, gene_name, donor_key, donor_str,
        n_donors, n_outcome_genes, run_hash, analysis_hash, out_file, out_Mfile
    )],
    file = file.path(out_dir, "mediation_run_index.tsv"),
    sep = "\t"
)
data.table::fwrite(
    run_info_dt[, .(
        iter, run_label, med_gene_id, gene_id, gene_name, donor_key, donor_str,
        n_donors, n_outcome_genes, run_hash, analysis_hash, out_file, out_Mfile
    )],
    file = file.path(out_dir, sprintf("mediation_run_index.%s.tsv", run_hash)),
    sep = "\t"
)

summary_dt <- run_info_dt[, .(
    iter, run_hash, mediation_run = run_label, donor_key, n_donors,
    baselineDEGs, mediationCarrierDEGs, mediationMDEGs, medSigBaselineDEGs,
    mediated_n, pass_both, loss, gain, loss_file, gain_file, mediated_file
)]
data.table::fwrite(
    summary_dt,
    file = file.path(out_dir, "mediation_vs_baseline_summary.tsv.gz"),
    sep = "\t",
    compress = "gzip"
)
data.table::fwrite(
    summary_dt,
    file = file.path(out_dir, sprintf("mediation_vs_baseline_summary.%s.tsv.gz", run_hash)),
    sep = "\t",
    compress = "gzip"
)

history_file <- file.path(out_dir, "run_history.tsv")
history_row <- data.table::data.table(
    run_hash = run_hash,
    run_started_at = run_started_at,
    run_finished_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    mediation = opt$mediation,
    FDRthr = FDRthr,
    medSelFDR = medSelFDR,
    forceSingleBaseline = forceSingleBaseline,
    n_mediators = nrow(mediator_plan_dt),
    n_mediators_eligible = nrow(mediator_plan_run),
    n_completed = nrow(run_info_dt),
    n_unique_donor_sets = nrow(donor_sets),
    outcome_donor_key = outcome_donor_info$donor_key,
    outcome_donor_str = outcome_donor_info$donor_str,
    outcome_gene_hash = outcome_gene_hash,
    med_degs_hash = med_degs_hash,
    med_expr_hash = med_expr_hash
)
history_dt <- if (file.exists(history_file)) {
    data.table::rbindlist(list(data.table::fread(history_file), history_row), fill = TRUE)
} else {
    history_row
}
data.table::fwrite(history_dt, file = history_file, sep = "\t")
message(Sys.time(), sprintf(" - Run metadata recorded in %s", history_file))
#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
