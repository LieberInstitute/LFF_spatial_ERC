## Louise Huuki-Myers & Bernie Mulvey, June 2025
## Run VoomLmFit on ERC clusters

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("here")
library("data.table")
library("getopt")
library("tidyverse")
library("sessioninfo")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype <- "sn_broad"

print(opt)

if(opt$datatype == "Visium"){
    pb_fn <- here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS")
    batch <- "Visium_slide"
} else if(opt$datatype == "sn_broad"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS")
    batch <- "exp_round"
}else if(opt$datatype == "sn_fine"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS")
    batch <- "exp_round"
} else {
    stop("non-valid datatype")
}

message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, basename(pb_fn)))

#### Set up dirs ####
data_dir <- here("processed-data", "12_voomLmFit", "03_Clusterwise_voomLmFit_interaction", sprintf("vlmf_%s", opt$datatype))
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
sce_pb <- readRDS(pb_fn)
dim(sce_pb)
table(sce_pb$registration_variable)

clusters <- levels(sce_pb$registration_variable)
names(clusters) <- clusters

message(Sys.time(), " - Loop voomlmFit by cluster")

lmf_summary <- map_dfr(clusters, function(clus){
    
    dge <- sce_pb[,sce_pb$registration_variable == clus]

    des <- model.matrix(~APOE_carrier_syn*Anc_Afr + Sex + Age + pseudo_expr_chrM_ratio, data = colData(dge))
    des <- as.data.frame(des)
    
    # filter low expression genes
    dge <- edgeR::calcNormFactors(dge)
    keep <- edgeR::filterByExpr.DGEList(dge,design=des)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- edgeR::calcNormFactors(dge)
    
    message(Sys.time(), sprintf(" - voomLmFit - cluster: %s, block= '%s', ncol: %s, ngene: %i", clus, batch, ncol(dge), nrow(dge$genes)))
    
    # make these more readable
    colnames(des) <- gsub(colnames(des),pattern="_syn",replacement="_")
    
    ## run voomLmFit for the pseudobulked data, referring donor to duplicateCorrelation; 
    ## using an adaptive span (number of genes, based on the number of genes in the dge) for smoothing the mean-variance trend
    v.swt <- voomLmFit(dge,design = des,block = as.factor(dge$samples[[batch]]),adaptive.span = T,sample.weights = T)
    
    v.swt.fit.e <- eBayes(v.swt)
    
    "APOE_carrier_E4:Anc_Afr"
    
    v.swt.e.tt <- topTable(v.swt.fit.e, 
                           coef = "APOE_carrier_E4:Anc_Afr",
                           number=Inf, 
                           adjust.method = "BH") |>
        mutate(data_type = opt$datatype, 
               cluster = clus, 
               .before = 1) |>
        arrange(adj.P.Val)
    
    message("Done - Save data")
    saveRDS(v.swt.e.tt, file = here(data_dir, sprintf("voomLmFit_%s_%s.rds", opt$datatype, clus)))
    return(purrr::map_int(v.swt.e.tt, ~sum(.x$adj.P.Val < 0.05)))
})

lmf_summary <- lmf_summary |>
    add_column(cluster = clusters, .before=1)

write.csv(lmf_summary, file = here(data_dir, sprintf("vlmf_FDR05_summary-%s.csv", opt$datatype)), row.names = FALSE)

# slurmjobs::job_single('03_Clusterwise_voomLmFit_interaction_sn_broad', create_shell = TRUE, memory = '25G', command = "Rscript 03_Clusterwise_voomLmFit_interaction.R --datatype sn_broad")
# slurmjobs::job_single('03_Clusterwise_voomLmFit_interaction_sn_fine', create_shell = TRUE, memory = '10G', command = "Rscript 03_Clusterwise_voomLmFit_interaction.R --datatype sn_fine")
# slurmjobs::job_single('03_Clusterwise_voomLmFit_interaction_Visium', create_shell = TRUE, memory = '25G', command = "Rscript 03_Clusterwise_voomLmFit_interaction.R --datatype Visium")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()