## Louise Huuki-Myers, April 2025
## Examine covaraites in Visium dataset for DGE
## adapted from https://github.com/LieberInstitute/dlpfc_asd/blob/a250a1d7e20bd754c5f1186aa96ce0752d55e556/code/08_pseudoBulkDGE_s/02_covariate_analysis.R

#### Set Up ####
library("spatialLIBD")
library("SingleCellExperiment")
# library("scran")
# library("BayesSpace")
library("tidyverse")
library("ggpubr")
library("ggrepel")
library("here")
library("sessioninfo")
library("variancePartition")
library("getopt")

data_dir <- here("plots", "09_pseudoBulkDGE_Visium", "02_VariancePartition_Visium")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "09_pseudoBulkDGE_Visium", "02_VariancePartition_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
scec <- matrix(
    c("spd", "s", "1", "character", "Name of SpD to test"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

# opt <- list(spd = "Vasc_Sp09D08")

#### Load the data ####
spe_pb <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))
stopifnot(opt$spd %in% levels(spe_pb$SpD_syn))

## subset 
spe_pb <- spe_pb[,spe_pb$SpD_syn == opt$spd]
message(sprintf("Subset to %s - %i spots", opt$spd, nrow(spe_pb)))

# Extract pheontype data
pd <- as.data.frame(colData(spe_pb)) 
           
#### Assess correlation between variables ####
form <- ~ APOE + APOE_num + sample_id + Sex + Age + Anc_Afr + Rin + Visium_slide + round + ncells + pseudo_sum_umi + pseudo_expr_chrM + pseudo_expr_chrM_ratio

C <- canCorPairs(form, pd)

pdf(here(plot_dir, sprintf("Visium_variable_cor_matrix-%s.pdf", opt$spd)))
plotCorrMatrix(C)
title(opt$spd)
dev.off()

#### Variance Partition ####

my_forms <- list(
    global = ~ (1 | APOE_carrier) + APOE_num + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium_slide) + (1 | round) + ncells,
    carrier_i = ~ (1 | APOE_carrier:Anc_Afr) + (1 | APOE_carrier:Sex) + (1 | APOE_carrier) + (1 | APOE_carrier) + Anc_Afr + Age + Rin + (1 | Visium_slide) + (1 | round) + ncells,
    numeric_i = ~ APOE_num:Anc_Afr + (1 | APOE_num:Sex) + APOE_num + (1 | Sex) + Age + Rin + (1 | Visium_slide) + (1 | round) + ncells
)

varPart_summary <- map2_dfr(my_forms, names(my_forms), function(form, form_name){
    
    message(Sys.time(), " - VarPart: ", form_name)
    varPart <- fitExtractVarPartModel(logcounts(spe_pb), form, pd)
    
    varPart_summary <- as.data.frame(apply(varPart, 2, summary)) |>
        rownames_to_column("metric") |>
        pivot_longer(!metric) |>
        mutate(form = form_name)
    
    vp <- sortCols(varPart)
    
    ## violin plot 
    sn_vp_violin <- plotVarPart(vp) + labs(title = opt$spd, subtitle = form_name)
    ggsave(sn_vp_violin, filename = here(plot_dir, sprintf("Visium_VarPart_violin_%s_-%s.png", form_name, opt$spd)))
    
    return(varPart_summary)
})

write_csv(varPart_summary, file = here(data_dir, sprintf("Visium_varPart_summary_%s.csv", opt$spd)))

# slurmjobs::job_single('02_VariancePartition_Visium', create_shell = TRUE, memory = '25G', command = "Rscript 02_VariancePartition_Visium.R")

# slurmjobs::job_loop(loops = list(spd = levels(spe_pb$SpD_syn)),
#                     name = "02_VariancePartition_Visium",
#                     create_shell = TRUE,
#                     create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()