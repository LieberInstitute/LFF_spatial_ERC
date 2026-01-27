#!/bin/env Rscript

## Louise Huuki-Myers & Bernie Mulvey, June 2025
## ERC-only mediation screen (Astro mediators -> Oligo.3 outcomes)

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("SummarizedExperiment")
library("BiocParallel")
library("here")
library("getopt")
library("sessioninfo")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- tryCatch(getopt(scec), error = function(e) list())
if (is.null(opt$datatype) || is.na(opt$datatype)) opt$datatype <- "sn_fine"

if (opt$datatype != "sn_fine") {
    stop("This mediation script expects sn_fine pseudo-bulk input (Astro.1/2/3 + Oligo.3).")
}

workers <- as.integer(Sys.getenv("MEDIATION_WORKERS", "10"))
if (is.na(workers) || workers < 1) workers <- 1

message(Sys.time(), sprintf(" - Datatype = %s (workers = %i)", opt$datatype, workers))

#### Paths ####
ddir <- here("processed-data")
pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn", "sce_pseudo_DGE-cell_type_anno.RDS")
degs_fn <- file.path(ddir, "13_compile_DGE/01_compile_DGE/sn_fine", "DGE_results_carrier_sn_fine.Rds")

out1_dir <- here("code", "22_Mediation", "out1")
if (!dir.exists(out1_dir)) dir.create(out1_dir, recursive = TRUE)

#### Load data ####
sce_pb <- readRDS(pb_fn)
sce_pb <- sce_pb[, sce_pb$BrNum != "Br1289"]

pdge <- readRDS(degs_fn)
degs <- as.data.table(pdge)[
    grepl("^Astro", cluster) & vlmf_adj.P.Val < 0.05 & contrast == "carrier",
    .(cluster, gene_id, gene_name, gene_type, vlmf_logFC, vlmf_AveExpr, vlmf_t, vlmf_P.Value, vlmf_adj.P.Val, cell_type_broad)
]
setorder(degs, cluster, gene_id)
degs[, iter := .I]

message(Sys.time(), sprintf(" - Astro DEGs for mediation: %i", nrow(degs)))

#### Outcome cluster setup ####
out_clus <- "Oligo.3"
batch <- "exp_round"

dge_base <- sce_pb[, sce_pb$registration_variable == out_clus]

baseline_des <- model.matrix(
    #~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,

    data = colData(dge_base)
)
baseline_des <- as.data.frame(baseline_des)

## Filter low expression genes once (baseline design)
dge_base <- edgeR::calcNormFactors(dge_base)
keep <- edgeR::filterByExpr(dge_base, design = baseline_des)
dge_base <- dge_base[keep, , keep.lib.sizes = FALSE]
dge_base <- edgeR::calcNormFactors(dge_base)

run_one <- function(i) {
    med_row <- degs[i]
    med_clus <- as.character(med_row$cluster)
    gene_id <- as.character(med_row$gene_id)
    gene_name <- as.character(med_row$gene_name)
    run_label <- sprintf("%s|%s|%s", med_clus, gene_id, gene_name)

    file_label <- gsub("\\|", "_", run_label)
    out_file <- file.path(out1_dir, sprintf("%s_%02d.tab.gz", file_label, i))

    # Get mediator expression
    med_sce <- sce_pb[gene_id, as.character(sce_pb$registration_variable) == med_clus]
    med_brnum <- as.character(colData(med_sce)$BrNum)

    # Find common donors between outcome (dge_base) and mediator cluster
    out_brnum <- as.character(dge_base$samples$BrNum)
    common_donors <- intersect(out_brnum, med_brnum)

    if (length(common_donors) < 10) {
        warning(sprintf("Too few common donors (%d) for %s", length(common_donors), run_label))
        return(NULL)
    }

    # Subset outcome data to common donors
    dge_sub <- dge_base[, out_brnum %in% common_donors]

    # Get mediator values aligned to outcome donors
    med_vec_full <- as.numeric(assay(med_sce, "logcounts")[1, ])
    names(med_vec_full) <- med_brnum
    med_vec <- med_vec_full[as.character(dge_sub$samples$BrNum)]

    if (anyNA(med_vec)) {
        stop(sprintf("Missing mediator values for %s (%s)", gene_id, med_clus))
    }

    des <- model.matrix(
        ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio + med_vec,
        data = data.frame(dge_sub$samples, med_vec = med_vec)
    )
    des <- as.data.frame(des)
    colnames(des) <- gsub(colnames(des), pattern = "_syn", replacement = "_")

    v.swt <- voomLmFit(
        dge_sub,
        design = des,
        block = as.factor(dge_sub$samples[[batch]]),
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
    tt[, `:=`(run = run_label, data_type = opt$datatype, cluster = out_clus, contrast = "carrier",
              n_donors = length(common_donors))]
    setcolorder(tt, c("run", "data_type", "cluster", "contrast", "outcome_gene_id", "n_donors",
                      setdiff(names(tt), c("run", "data_type", "cluster", "contrast", "outcome_gene_id", "n_donors"))))

    data.table::fwrite(tt, file = out_file, sep = "\t", compress = "gzip")
    out_file
}

message(Sys.time(), " - Loop1: cluster-specific astro mediators")
bpparam <- BiocParallel::MulticoreParam(workers = workers)
out_files <- BiocParallel::bplapply(seq_len(nrow(degs)), run_one, BPPARAM = bpparam)
message(Sys.time(), sprintf(" - Loop1 complete (%i files)", length(out_files)))

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
