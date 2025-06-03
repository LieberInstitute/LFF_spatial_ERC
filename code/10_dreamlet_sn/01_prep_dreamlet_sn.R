## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "01_prep_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

colnames(sce) <- paste0(sce$sample_id, "_", sce$Barcode)
rownames(sce) <- rowData(sce)$gene_name

## duplicated gene_names
sce <- sce[!duplicated(rownames(sce)),]

AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    dplyr::filter(symbol %in% rowData(sce)$gene_name) 

#### Aggregate to pseudobulk ####
message(Sys.time(), ' - Aggregate to pseudobulk')
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",
                            cluster_id = "cell_type_broad",
                            sample_id = "sample_id",
                            verbose = FALSE
)

## one assay per cell type
assayNames(pb)

dim(pb)

# Normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, ~APOE, min.count = 5)

# the resulting object of class dreamletProcessedData stores
# normalized data and other information
res.proc

details(res.proc)

# view details of dropping samples
details(res.proc)

## save 
saveRDS(res.proc, file = here(data_dir, "sn_res_proc.rds"))

# show voom plot for each cell clusters
pdf(here(plot_dir, "sn_dreamlet_voom.pdf"))
plotVoom(res.proc)
dev.off()

##### run variance partitioning ####
message(Sys.time(), ' - variance partition')

# my_forms <- list(
#     carrier = ~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
#     apoe = ~ (1 | APOE) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
#     e4e4 = ~ (1 | APOE_E4E4) + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | exp_round)  + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
#     # interaction syntax x + (x | g)
#     carrier_i = ~ Anc_Afr + (Anc_Afr | APOE_carrier) + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
#     apoe_i = ~ Anc_Afr + (Anc_Afr | APOE) + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
#     e4e4_i = ~ Anc_Afr + (Anc_Afr | APOE_E4E4) + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells
# )

vp.lst <- fitVarPart(res.proc,~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round)) #subsets_Mito_percent + ncells

saveRDS(vp.lst, file = here(data_dir, "dreamlet_fitVarPart_sn.RDS"))
# vp.lst <- readRDS(here(data_dir, "dreamlet_fitVarPart_sn.RDS"))

# Show variance fractions at the gene-level for each cell type
# genes <- vp.lst$gene[2:4]

vp.lst[order(vp.lst$APOE_carrier, decreasing = TRUE)[1:100],]

top20_APOE_varPart <- vp.lst$gene[order(vp.lst$APOE_carrier, decreasing = TRUE)[1:20]]

pdf(here(plot_dir, "sn_dreamlet_VarPart.pdf"), height = 11)
plotVarPart(vp.lst, label.angle = 45)
plotPercentBars(vp.lst[vp.lst$gene %in% AD_risk$symbol, ]) + labs(title = "Risk Genes")
plotPercentBars(vp.lst[vp.lst$gene %in% top20_APOE_varPart, ]) + labs(title = "Top20 APOE")
dev.off()

# slurmjobs::job_single('01_prep_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 01_prep_dreamlet_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

