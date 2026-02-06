#!/bin/env Rscript

## Louise Huuki-Myers & Bernie Mulvey, June 2025
## adapting for mediation screen by Geo Pertea, January 2026
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

# Helpers
md5_string <- function(x) {
    tf <- tempfile()
    writeLines(x, con = tf, useBytes = TRUE)
    md5 <- as.character(tools::md5sum(tf))
    unlink(tf)
    md5
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
    lines <- c(
        sprintf("run: %s", info$run_label),
        sprintf("med_gene_id: %s", info$med_gene_id),
        sprintf("gene_id: %s", info$gene_id),
        sprintf("gene_name: %s", info$gene_name),
        sprintf("donor_key: %s", info$donor_key),
        sprintf("n_donors: %d", info$n_donors),
        sprintf("n_outcome_genes: %d", info$n_outcome_genes),
        sprintf("design_formula: %s", info$design_formula),
        sprintf("design_dim: %dx%d", info$design_dim[1], info$design_dim[2]),
        sprintf("donors: %s", info$donor_str)
    )
    writeLines(lines, con = file)
}

read_run_info <- function(file) {
    if (!file.exists(file)) return(NULL)
    lines <- readLines(file, warn = FALSE)
    kv <- strsplit(lines, ": ", fixed = TRUE)
    keys <- vapply(kv, `[[`, character(1), 1)
    vals <- vapply(kv, function(x) paste(x[-1], collapse = ": "), character(1))
    donor_str <- vals[keys == "donors"]
    donors <- if (length(donor_str) > 0 && nzchar(donor_str)) strsplit(donor_str, ",", fixed = TRUE)[[1]] else character(0)
    list(
        run_label = vals[keys == "run"],
        med_gene_id = vals[keys == "med_gene_id"],
        gene_id = vals[keys == "gene_id"],
        gene_name = vals[keys == "gene_name"],
        donor_key = vals[keys == "donor_key"],
        donor_str = donor_str,
        donors = donors,
        n_donors = as.integer(vals[keys == "n_donors"]),
        n_outcome_genes = as.integer(vals[keys == "n_outcome_genes"])
    )
}

write_fdr_table <- function(dt, file) {
    ## write empty gzip with header when no rows to avoid corrupt gzip output
    if (nrow(dt) < 1) {
        con <- gzfile(file, "wt")
        on.exit(close(con), add = TRUE)
        if (ncol(dt) > 0) {
            writeLines(paste(names(dt), collapse = "\t"), con = con, useBytes = TRUE)
        }
        return(invisible(NULL))
    }
    data.table::fwrite(dt, file = file, sep = "\t", compress = "gzip")
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
workers = 10

test_rows <- integer(0)
#test_rows <- c(1,10,20) ## change and uncomment to test specific mediator rows
if (length(test_rows) > 0) {
        message("Mediation test rows provided: ", paste(test_rows, collapse = ", "))
}

message(Sys.time(), sprintf(" - Mediation type = %s (workers = %i)", opt$mediation, workers))

#### Paths ####
ddir <- here("processed-data")
mdir <- here("processed-data", "22_Mediation")

pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn", "sce_pseudo_DGE-cell_type_anno.RDS")

out_dir <- here("code", "22_Mediation", sprintf("out-%s", opt$mediation))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
baseline_dir <- file.path(out_dir, "baseline")
dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)

#### Load X->Y pseudo bulk data ####
sce_pb <- readRDS(pb_fn)
## always drop Br1289
sce_pb <- sce_pb[, sce_pb$BrNum != "Br1289"]

#### Outcome cluster setup (same for every scenario) ####
## we can use dge_base colData() for all donor metadata

out_clus <- "Oligo.3"
batch <- "exp_round"
sce_base <- sce_pb[, sce_pb$registration_variable == out_clus]
donor_meta <- as.data.frame(colData(sce_base))
baseline_des <- model.matrix(
    baseline_formula,
    data = donor_meta
)
baseline_des <- as.data.frame(baseline_des)


med_degs <- NULL ## placeholder for the selected mediator DEGs, gene metadata
med_expr <- NULL ## placeholder for the mediator gene expression data (logcounts)
## these should be loaded depending on mediation type

if (opt$mediation == "erc_astro") { ## "green" design 3: Mediator = ERC Astro DEGs
  degs_fn <- file.path(ddir, "13_compile_DGE/01_compile_DGE/sn_fine", "DGE_results_carrier_sn_fine.Rds")
  pdge <- readRDS(degs_fn)
  med_degs <- as.data.table(pdge)[
      grepl("^Astro", cluster) & vlmf_adj.P.Val < 0.05 & contrast == "carrier",
      .(cluster, gene_id, gene_name, gene_type, vlmf_logFC, vlmf_AveExpr, vlmf_t, vlmf_P.Value, vlmf_adj.P.Val, cell_type_broad)
  ]
  setorder(med_degs, cluster, gene_id)
  med_degs[, med_gene_id := paste0(cluster, "|", gene_id)] ## unique mediator ids, matches med_expr rownames
  med_degs[, run := sprintf("%s|%s|%s", cluster, gene_id, gene_name)]
  med_degs[, iter := .I]

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
  message(Sys.time(), sprintf(" - ERC Astro DEGs for mediation: %i", nrow(med_degs)))
  # fwrite med_degs and med_expr for "erc_astro" mediation in out_dir if they do not exist
  fn_med_degs <- file.path(out_dir, "med_degs.qs2")
  fn_med_expr <- file.path(out_dir, "med_expr.qs2")
  if (!file.exists(fn_med_degs)) {
      qs_save(med_degs, fn_med_degs)
  }
  if (!file.exists(fn_med_expr)) {
      qs_save(med_expr, fn_med_expr)
  }
} else if (opt$mediation == "lc_astro") { ## "blue" design 2: Mediator = LC Astro DEGs from LC spatial domains
    ## Load LC Astro DEGs (FDR < 0.05)
    f_lca_dge <- file.path(mdir, "LC_Astro_E4vE2_DGEout_Bernie.rds")
    if (!file.exists(f_lca_dge)) {
        stop(sprintf("LC Astro DEG file not found: %s", f_lca_dge))
    }
    lca_dge <- as.data.table(readRDS(f_lca_dge))
    med_degs <- lca_dge[adj.P.Val < 0.05]
    if (nrow(med_degs) < 1) {
        stop("No LC Astro DEGs with FDR < 0.05 found.")
    }
    ## Standardize column names and add required fields
    med_degs[, cluster := "LCAstro"]
    med_degs[, med_gene_id := paste0("LCAstro|", gene_id)]
    med_degs[, run := sprintf("LCAstro|%s|%s", gene_id, gene_name)]
    med_degs[, iter := .I]
    setorder(med_degs, gene_id)
    message(Sys.time(), sprintf(" - LC Astro DEGs for mediation: %i", nrow(med_degs)))

    ## Load LC Astro pseudobulk counts (DGEList)
    f_lca_pbcounts <- file.path(mdir, "lc_astro_pseudobulk_counts_by_donor.rds")
    if (!file.exists(f_lca_pbcounts)) {
        stop(sprintf("LC Astro pseudobulk counts file not found: %s", f_lca_pbcounts))
    }
    dgl_lca <- readRDS(f_lca_pbcounts)
    ## Determine overlapping donors between LC Astro and ERC Oligo.3 outcome
    ## after leave-one-out testing, let's drop BR5529 from dgl_lca
    dgl_lca <- dgl_lca[, !colnames(dgl_lca) %in% "Br5529"] ## or also "Br6423"
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

    ## Compute log-CPM
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

    ## Cache med_degs and med_expr
    fn_med_degs <- file.path(out_dir, "med_degs.qs2")
    fn_med_expr <- file.path(out_dir, "med_expr.qs2")
    if (!file.exists(fn_med_degs)) {
        qs_save(med_degs, fn_med_degs)
    }
    if (!file.exists(fn_med_expr)) {
        qs_save(med_expr, fn_med_expr)
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

## subset dge_base to the column names in med_expr -- the actual donors to keep
dge_base <- sce_base[, sce_base$BrNum %in% colnames(med_expr)]

## Filter low expression genes once (baseline design)
fn_dge_base <- file.path(out_dir, "dge_base.qs2")
dge_base <- edgeR::calcNormFactors(dge_base)
keep <- edgeR::filterByExpr(dge_base, design = baseline_des)
dge_base <- dge_base[keep, , keep.lib.sizes = FALSE]
dge_base <- edgeR::calcNormFactors(dge_base)
## save the file as qs2
qs_save(dge_base, fn_dge_base)

run_one <- function(i) {
    med_row <- med_degs[i] ## mediator row to use as covariate
    gene_id <- as.character(med_row$gene_id)
    gene_name <- as.character(med_row$gene_name)
    run_label <- as.character(med_row$run)
    med_gene_id <- as.character(med_row$med_gene_id)
    file_label <- gsub("\\|", "_", run_label)
    out_file <- file.path(out_dir, sprintf("%02d_%s.tab.gz", i, file_label)) ## topTable for carrier contrast
    out_Mfile <- file.path(out_dir, sprintf("%02d_%s.M.tab.gz", i, file_label)) ## topTable for med_vec contrast
    sufdr <- stringr::str_split_1(as.character(FDRthr), "\\.")[[2]]
    out_fdr_file <- sub("\\.tab\\.gz$", sprintf(".fdr%s.tab.gz", sufdr), out_file) ## carrier contrast FDR<0.05
    out_Mfdr_file <- sub("\\.tab\\.gz$", sprintf(".fdr%s.tab.gz", sufdr), out_Mfile) ## med_vec contrast FDR<0.05
    run_info_file <- file.path(out_dir, sprintf("%02d_%s.txt", i, file_label))
    if (file.exists(out_file) && file.exists(out_Mfile) && file.exists(run_info_file)) {
        cached <- read_run_info(run_info_file)
        if (!is.null(cached)) {
            ## backfill missing fdr outputs for cached runs
            if (!file.exists(out_fdr_file)) {
                tt_cached <- data.table::fread(out_file)
                tt_fdr <- tt_cached[adj.P.Val < FDRthr]
                write_fdr_table(tt_fdr, out_fdr_file)
            }
            if (!file.exists(out_Mfdr_file)) {
                ttM_cached <- data.table::fread(out_Mfile)
                ttM_fdr <- ttM_cached[adj.P.Val < FDRthr]
                write_fdr_table(ttM_fdr, out_Mfdr_file)
            }
            message(Sys.time(), sprintf(" - Using cached run: %s", out_file))
            return(list(
                out_file = out_file,
                out_Mfile = out_Mfile,
                run_label = cached$run_label,
                med_gene_id = cached$med_gene_id,
                gene_id = cached$gene_id,
                gene_name = cached$gene_name,
                donor_key = cached$donor_key,
                donor_str = cached$donor_str,
                donors = list(cached$donors),
                n_donors = cached$n_donors,
                n_outcome_genes = cached$n_outcome_genes
            ))
        }
    }
    ## Get mediator expression from precomputed matrix, drop NA donors
    if (!(med_gene_id %in% rownames(med_expr))) {
        stop(sprintf("med_gene_id not found in med_expr: %s", med_gene_id))
    }
    med_vec_full <- as.numeric(med_expr[med_gene_id, ])
    names(med_vec_full) <- colnames(med_expr)
    med_vec_full <- med_vec_full[!is.na(med_vec_full)]
    med_brnum <- names(med_vec_full)

    ## Important: find common donors between outcome (dge_base) and mediator cluster
    out_brnum <- as.character(dge_base$samples$BrNum)
    common_donors <- intersect(out_brnum, med_brnum)

    if (length(common_donors) < 10) {
        warning(sprintf("Too few common donors (%d) for %s", length(common_donors), run_label))
        return(NULL)
    }

    donor_info <- make_donor_key(common_donors)

    # Subset outcome data to common donors
    dge_sub <- dge_base[, out_brnum %in% common_donors]

    # Get mediator values aligned to outcome donors
    med_vec <- med_vec_full[as.character(dge_sub$samples$BrNum)]

    if (anyNA(med_vec)) {
        stop(sprintf("Missing mediator values for %s (%s)", gene_id, run_label))
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

    ## add identity contrast for the med_vec coefficient alongside carrier
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
              n_donors = length(common_donors))]
    setcolorder(tt, c("run", "cluster", "contrast", "outcome_gene_id", "n_donors",
                      setdiff(names(tt), c("run", "cluster", "contrast", "outcome_gene_id", "n_donors"))))
    ## mediator contrast table
    ttM <- topTable(v.swt.fit.e, coef = "med_vec", number = Inf, adjust.method = "BH")
    ttM <- as.data.table(ttM, keep.rownames = "outcome_gene_id")
    ttM[, `:=`(run = run_label, cluster = out_clus, contrast = "med_vec",
               n_donors = length(common_donors))]
    setcolorder(ttM, c("run", "cluster", "contrast", "outcome_gene_id", "n_donors",
                       setdiff(names(ttM), c("run", "cluster", "contrast", "outcome_gene_id", "n_donors"))))

    ## write carrier and mediator results plus fdr-filtered subsets
    data.table::fwrite(tt, file = out_file, sep = "\t", compress = "gzip")
    write_fdr_table(tt[adj.P.Val < FDRthr], out_fdr_file)
    data.table::fwrite(ttM, file = out_Mfile, sep = "\t", compress = "gzip")
    write_fdr_table(ttM[adj.P.Val < FDRthr], out_Mfdr_file)

    write_run_info(run_info_file, list(
        run_label = run_label,
        med_gene_id = med_gene_id,
        gene_id = gene_id,
        gene_name = gene_name,
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str,
        n_donors = length(common_donors),
        n_outcome_genes = nrow(dge_sub),
        design_formula = mediation_formula_str,
        design_dim = c(nrow(des), ncol(des))
    ))

    list(
        out_file = out_file,
        out_Mfile = out_Mfile,
        run_label = run_label,
        med_gene_id = med_gene_id,
        gene_id = gene_id,
        gene_name = gene_name,
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str,
        donors = list(donor_info$donors),
        n_donors = length(common_donors),
        n_outcome_genes = nrow(dge_sub)
    )
}

message(Sys.time(), sprintf(" - Starting %d workers for %d %s mediation runs",
                                               workers, nrow(med_degs), opt$mediation))
bpparam <- BiocParallel::MulticoreParam(workers = workers)
run_results <- BiocParallel::bplapply(seq_len(nrow(med_degs)), run_one, BPPARAM = bpparam)
run_results <- Filter(Negate(is.null), run_results)
message(Sys.time(), sprintf(" - Runs completed (%i files)", length(run_results)))

run_info_dt <- data.table::rbindlist(run_results, fill = TRUE)

if (nrow(run_info_dt) < 1) {
    stop("No mediation runs completed successfully.")
}

## Unique donor sets
donor_sets <- unique(run_info_dt[, .(donor_key, donor_str, donors)], by = "donor_key")

if (forceSingleBaseline) {
    # Find the donor set with maximum donors
    max_idx <- which.max(sapply(donor_sets$donors, length))
    donor_sets_to_run <- donor_sets[max_idx]
    single_baseline_key <- donor_sets_to_run$donor_key
    message(Sys.time(), sprintf(" - Using single baseline with %d donors (forceSingleBaseline=TRUE)",
                                 length(donor_sets_to_run$donors[[1]])))
} else {
    donor_sets_to_run <- donor_sets
    single_baseline_key <- NULL
}

run_baseline <- function(donor_key, donors, donor_str) {
    out_file <- file.path(baseline_dir, sprintf("baseline_%s.tab.gz", donor_key))
    sufdr <- stringr::str_split_1(as.character(FDRthr), "\\.")[[2]]
    out_fdr_file <- file.path(baseline_dir, sprintf("baseline_%s.fdr%s.tab.gz", donor_key, sufdr))
    if (file.exists(out_file)) {
        message(Sys.time(), sprintf(" - Using cached baseline: %s", out_file))
        tt <- data.table::fread(out_file)
        baseline_sig <- tt[adj.P.Val < FDRthr, outcome_gene_id]
        if (!file.exists(out_fdr_file)) {
            data.table::fwrite(tt[adj.P.Val < FDRthr], file = out_fdr_file, sep = "\t", compress = "gzip")
        }
        return(list(
            donor_key = donor_key,
            donor_str = donor_str,
            n_donors = length(donors),
            out_file = out_file,
            sig_genes = list(baseline_sig)
        ))
    }

    out_brnum <- as.character(dge_base$samples$BrNum)
    common_donors <- intersect(out_brnum, donors)
    if (length(common_donors) < 10) {
        warning(sprintf("Too few common donors (%d) for baseline %s", length(common_donors), donor_key))
        return(NULL)
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
    data.table::fwrite(tt[adj.P.Val < FDRthr], file = out_fdr_file, sep = "\t", compress = "gzip")

    baseline_sig <- tt[adj.P.Val < FDRthr, outcome_gene_id]
    list(
        donor_key = donor_key,
        donor_str = donor_str,
        n_donors = length(common_donors),
        out_file = out_file,
        sig_genes = list(baseline_sig)
    )
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

# When forceSingleBaseline=TRUE, map all donor_keys to the single baseline
if (forceSingleBaseline && !is.null(single_baseline_key)) {
    single_baseline_genes <- baseline_map[[single_baseline_key]]
    for (dk in donor_sets$donor_key) {
        baseline_map[[dk]] <- single_baseline_genes
    }
}

data.table::fwrite(
    baseline_dt[, .(donor_key, donor_str, n_donors, out_file)],
    file = file.path(out_dir, "baseline_index.tsv"),
    sep = "\t"
)

gene_name_vec <- NULL
if (inherits(dge_base, "DGEList")) {
    if (!is.null(dge_base$genes) && "gene_name" %in% colnames(dge_base$genes)) {
        gene_name_vec <- as.character(dge_base$genes$gene_name)
    }
} else if ("gene_name" %in% colnames(rowData(dge_base))) {
    gene_name_vec <- as.character(rowData(dge_base)$gene_name)
}
if (is.null(gene_name_vec)) {
    gene_name_vec <- rep(NA_character_, nrow(dge_base))
}
out_gene_map <- data.table::data.table(
    outcome_gene_id = rownames(dge_base),
    outcome_gene_name = gene_name_vec
)
out_gene_map <- unique(out_gene_map, by = "outcome_gene_id")

summary_list <- lapply(seq_len(nrow(run_info_dt)), function(i) {
    run_row <- run_info_dt[i]
    ## read carrier and mediator results for step 3 evaluation
    tt <- data.table::fread(run_row$out_file)
    ttM <- data.table::fread(run_row$out_Mfile)
    carrier_sig <- tt[adj.P.Val < FDRthr, outcome_gene_id]
    med_sig <- ttM[adj.P.Val < FDRthr, outcome_gene_id]
    base_genes <- baseline_map[[run_row$donor_key]]
    if (is.null(base_genes)) {
        stop(sprintf("Missing baseline for donor_key %s", run_row$donor_key))
    }
    ## step 3b attenuation and step 3a mediator significance (baseline-only)
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

    file_label <- gsub("\\|", "_", run_row$run_label)
    loss_file <- file.path(out_dir, sprintf("%s.loss.txt", file_label))
    gain_file <- file.path(out_dir, sprintf("%s.gain.txt", file_label))
    mediated_file <- file.path(out_dir, sprintf("%s.mediated.txt", file_label))
    lost_labels <- sort(as.character(lost_labels))
    gained_labels <- sort(as.character(gained_labels))
    mediated_labels <- sort(as.character(mediated_labels))
    lost_labels <- lost_labels[!is.na(lost_labels)]
    gained_labels <- gained_labels[!is.na(gained_labels)]
    mediated_labels <- mediated_labels[!is.na(mediated_labels)]
    ## write per-run gene lists
    writeLines(lost_labels, con = loss_file)
    writeLines(gained_labels, con = gain_file)
    writeLines(mediated_labels, con = mediated_file)

    data.table::data.table(
        mediation_run = run_row$run_label,
        donor_key = run_row$donor_key,
        n_donors = run_row$n_donors,
        baselineDEGs = length(base_genes),
        mediationCarrierDEGs = length(carrier_sig), ## carrier-significant outcome genes after adding mediator
        mediationMDEGs = length(med_sig), ## mediator-significant outcome genes in joint model
        medSigBaselineDEGs = length(step3a_genes), ## mediator-significant genes restricted to baseline DEG set
        mediated_n = length(mediated),
        pass_both = length(mediated) > 0,
        loss = length(lost),
        gain = length(gained),
        loss_file = loss_file,
        gain_file = gain_file,
        mediated_file = mediated_file
    )
})

summary_dt <- data.table::rbindlist(summary_list, fill = TRUE)
data.table::fwrite(
    summary_dt,
    file = file.path(out_dir, "mediation_vs_baseline_summary.tsv.gz"),
    sep = "\t",
    compress = "gzip"
)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
