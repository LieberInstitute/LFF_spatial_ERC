## Louise Huuki-Myers, May 2025
## Run Variance Partition on each cell types over several formulas

#### Set Up ####
library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")
library("variancePartition")
library("getopt")

data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "02_VariancePartition_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "08_pseudoBulkDGE_sn", "02_VariancePartition_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
scec <- matrix(
    c("cluster", "c", "1", "character", "Name of cluster",
      "cell_type", "ct", "1", "character", "Name of SpD to test"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

# opt <- list(cluster = "cell_type_broad", cell_type = "Excit")

#### Load the data ####
sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn", sprintf("sce_pseudo_DGE-%s.RDS", opt$cluster)))
stopifnot(opt$cell_type %in% levels(sce_pb[[opt$cluster]]))

## subset 
sce_pb <- sce_pb[,sce_pb[[opt$cluster]] == opt$cell_type]
message(sprintf("Subset to %s - %i samples", opt$cell_type, ncol(sce_pb)))

## add E4/E4 variable
sce_pb$APOE_E4E4 <- sce_pb$APOE == "E4/E4"
table(sce_pb$APOE_E4E4)

## scale sum_umi & Mitocondrial ratio
sce_pb$pseudo_sum_umi <- scale(sce_pb$pseudo_sum_umi)
sce_pb$pseudo_expr_chrM_ratio <- scale(sce_pb$pseudo_expr_chrM_ratio)

#### Assess correlation between variables ####
# Extract pheontype data
message(Sys.time(), " - canCorPairs")

pd <- as.data.frame(colData(sce_pb)) 

form <- ~ APOE + APOE_carrier + Sex + Age + Anc_Afr + Rin + exp_round + seq_round + pseudo_sum_umi + pseudo_expr_chrM_ratio + ncells

C <- canCorPairs(form, pd)

pdf(here(plot_dir, sprintf("sn_variable_cor_matrix-%s_%s.pdf", opt$cell_cluster, opt$cell_type)))
plotCorrMatrix(C)
title(opt$cell_type)
dev.off()

#### Variance Partition ####
my_forms <- list(
    carrier = ~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
    apoe = ~ (1 | APOE) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
    e4e4 = ~ (1 | APOE_E4E4) + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | exp_round)  + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
    # interaction syntax x + (x | g)
    carrier_i = ~ Anc_Afr + (Anc_Afr | APOE_carrier) + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
    apoe_i = ~ Anc_Afr + (Anc_Afr | APOE) + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells,
    e4e4_i = ~ Anc_Afr + (Anc_Afr | APOE_E4E4) + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round) + pseudo_expr_chrM_ratio + ncells
)

varPart_summary <- map2_dfr(my_forms, names(my_forms), function(form, form_name){
    
    message(Sys.time(), " - VarPart: ", form_name)
    
    ## filter genes?
    gkeep <- edgeR::filterByExpr(sce_pb, design=form, group=sce_pb$APOE)
    message("filter genes:", nrow(sce_pb) - length(gkeep))
    
    varPart <- fitExtractVarPartModel(logcounts(sce_pb), form, pd)
    
    varPart_summary <- as.data.frame(apply(varPart, 2, summary)) |>
        rownames_to_column("metric") |>
        pivot_longer(!metric) |>
        mutate(form = form_name,
               cell_type = opt$cell_type)
    
    vp <- sortCols(varPart)
    
    ## violin plot 
    sn_vp_violin <- plotVarPart(vp) + labs(title = opt$cell_type, subtitle = form_name)
    ggsave(sn_vp_violin, filename = here(plot_dir, sprintf("sn_VarPart_violin-%s_%s_-%s.png", opt$cluster, form_name, opt$cell_type)))
    
    return(varPart_summary)
})

write_csv(varPart_summary, file = here(data_dir, sprintf("sn_varPart_summary-%s_%s.csv", opt$cluster, opt$cell_type)))

# slurmjobs::job_single('02_VariancePartition_sn', create_shell = TRUE, memory = '25G', command = "Rscript 02_VariancePartition_sn.R")
# slurmjobs::job_loop(loops = list(cell_type = levels(sce_pb$cell_type_broad)),
#                     name = "02_VariancePartition_sn",
#                     create_shell = TRUE,
#                     create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()