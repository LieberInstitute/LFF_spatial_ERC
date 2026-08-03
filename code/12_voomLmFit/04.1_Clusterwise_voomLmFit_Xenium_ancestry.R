## Louise Huuki-Myers & Bernie Mulvey, Aug 2026
## Modified: ancestry-stratified voomLmFit on Xenium ERC clusters (E2 vs E4 within AA and EA)
## Goal: validate snRNA-seq ancestry-specific APOE findings in Xenium

#### Set Up ####

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("here")
library("data.table")
library("getopt")
library("tidyverse")
library("sessioninfo")

opt <- list()
opt$datatype <- "Xenium"

pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno.RDS")
batch <- "chip"

message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, basename(pb_fn)))

#### Set up dirs ####
data_dir <- here("processed-data", "12_voomLmFit", "04_Clusterwise_voomLmFit_ancestry", sprintf("vlmf_%s", opt$datatype))
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
sce_pb <- readRDS(pb_fn)

dim(sce_pb)
table(sce_pb$registration_variable)

## quick look at available metadata - check this against the assumptions below before trusting the full run
print(colnames(colData(sce_pb)))

## Clean up APOE carrier label (drop '+')
sce_pb$APOE_carrier_syn <- gsub("\\+", "", sce_pb$APOE_carrier)
table(sce_pb$APOE_carrier_syn)

## Ancestry variable check
table(sce_pb$Ancestry)

clusters <- levels(sce_pb$registration_variable)
names(clusters) <- clusters

## Add carrier x ancestry grouping
sce_pb$carrier_Anc <- paste0(sce_pb$APOE_carrier_syn, "_", sce_pb$Ancestry)
table(sce_pb$carrier_Anc)

carrier_levels <- sort(unique(sce_pb$APOE_carrier_syn))
anc_levels <- sort(unique(sce_pb$Ancestry))

message(Sys.time(), " - Loop voomlmFit by cluster (ancestry-stratified)")

lmf_summary <- map_dfr(clusters, possibly(function(clus){

    dge <- sce_pb[,sce_pb$registration_variable == clus]

    # table(dge$carrier_Anc)
    des <- model.matrix(~0 + carrier_Anc + Age,
                         data = colData(dge)) ## no Mito ratio for Xenium

    des <- as.data.frame(des)

    # filter low expression genes
    dge <- edgeR::calcNormFactors(dge)
    keep <- edgeR::filterByExpr.DGEList(dge, design = des)
    dge <- dge[keep,,keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge)

    message(Sys.time(), sprintf(" - voomLmFit - cluster: %s, block= '%s', ncol: %s, ngene: %i", clus, batch, ncol(dge), nrow(dge$genes)))

    # make these more readable
    colnames(des) <- gsub(colnames(des), pattern = "_syn", replacement = "_")

    ## run voomLmFit for the pseudobulked data, referring donor to duplicateCorrelation;
    ## using an adaptive span (number of genes, based on the number of genes in the dge) for smoothing the mean-variance trend
    v.swt <- voomLmFit(dge, design = des,
                        block = as.factor(dge$samples[[batch]]),
                        adaptive.span = T,
                        sample.weights = T)

    ## build contrasts from whatever the actual carrier_Anc levels are, rather than
    ## hardcoding "E2"/"E4"/"AA"/"EA" strings - protects against label mismatches
    coef1 <- paste0("carrier_Anc", carrier_levels[1], "_", anc_levels)  # e.g. E2_AA, E2_EA
    coef2 <- paste0("carrier_Anc", carrier_levels[2], "_", anc_levels)  # e.g. E4_AA, E4_EA

    if (!all(c(coef1, coef2) %in% colnames(des))) {
        stop(sprintf("Cluster %s: expected design columns not found. Looked for: %s. Design has: %s",
                      clus, paste(c(coef1, coef2), collapse = ", "), paste(colnames(des), collapse = ", ")))
    }

    cont <- makeContrasts(
        ## compare carrier by ancestry
        carrier_1 = sprintf("-%s + %s", coef1[1], coef2[1]),
        carrier_2 = sprintf("-%s + %s", coef1[2], coef2[2]),
        levels = des
    )
    colnames(cont) <- c(paste0("carrier_", anc_levels[1]), paste0("carrier_", anc_levels[2]))

    v.swt.fit <- contrasts.fit(v.swt, contrasts = cont)
    v.swt.fit.e <- eBayes(v.swt.fit)

    ## run top table over contrasts
    v.swt.e.tt <- purrr::map(colnames(cont), ~topTable(v.swt.fit.e, coef = .x, number = Inf, adjust.method = "BH") |>
                                 mutate(data_type = opt$datatype,
                                        cluster = clus,
                                        contrast = .x,
                                        .before = 1) |>
                                 arrange(adj.P.Val))

    names(v.swt.e.tt) <- colnames(cont)

    map(v.swt.e.tt, head)

    message("Done - Save data")
    saveRDS(v.swt.e.tt, file = here(data_dir, sprintf("voomLmFit_ancestry_%s_%s.rds", opt$datatype, clus)))

    return(purrr::map_int(v.swt.e.tt, ~sum(.x$adj.P.Val < 0.05)))

}, otherwise = NA))

lmf_summary <- lmf_summary |>
    add_column(cluster = clusters, .before = 1)

write.csv(lmf_summary, file = here(data_dir, sprintf("vlmf_ancestry_FDR05_summary-%s.csv", opt$datatype)), row.names = FALSE)

# slurmjobs::job_single('04.1_Clusterwise_voomLmFit_Xenium_ancestry', create_shell = TRUE, memory = '25G', command = "Rscript 04.1_Clusterwise_voomLmFit_Xenium_ancestry.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
