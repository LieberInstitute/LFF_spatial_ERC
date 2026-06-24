## Louise Huuki-Myers June 2026
## Run voomLmFit on Xenium data

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("here")
library("data.table")
library("getopt")
library("tidyverse")
library("sessioninfo")
library("qs2")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

opt$datatype <- "Xenium_Oligo.3_Astro"

if(opt$datatype == "Xenium_cell_type"){
    pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno.RDS")
} else if(opt$datatype == "Xenium_cell_type_anno_SpX"){
    pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno_SpX.RDS")
}else if(opt$datatype == "Xenium_Oligo.3_Astro"){
    pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-Oligo.3_Astro.RDS")
} else {
    stop("non-valid datatype")
}

batch = "chip"
message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, pb_fn))

#### Set up dirs ####
data_dir <- here("processed-data", "12_voomLmFit", "06_Clusterwise_voomLmFit_Xenium", sprintf("vlmf_%s", opt$datatype))
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
sce_pb <- readRDS(pb_fn)

dim(sce_pb)
table(sce_pb$registration_variable)

sce_pb$APOE_carrier_syn <- gsub("\\+", "", sce_pb$APOE_carrier)
table(sce_pb$APOE_carrier_syn)

clusters <- unique(sce_pb$registration_variable)
names(clusters) <- clusters

message(Sys.time(), " - Loop voomlmFit by cluster")

lmf_summary <- map_dfr(clusters, possibly(function(clus){
    
    dge <- sce_pb[,sce_pb$registration_variable ==clus]
    
    c_cell_type <- unique(dge$cell_type_anno)
    c_SpX <- unique(dge$SpX)

    des <- model.matrix(~APOE_carrier_syn + Age + Anc_Afr , data = colData(dge)) ## no Mito ratio for Xenium

    des <- as.data.frame(des)
    
    # filter low expression genes
    dge <- edgeR::calcNormFactors(dge)
    keep <- edgeR::filterByExpr.DGEList(dge,design=des)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- edgeR::calcNormFactors(dge)
    
    message(Sys.time(), sprintf(" - voomLmFit - cluster: %s, block= '%s', ncol: %s, ngene: %i", clus, batch, ncol(dge), nrow(dge$genes)))
    
    # make these more readable
    colnames(des) <- gsub(colnames(des), pattern="_syn", replacement="_")
    
    ## run voomLmFit for the pseudobulked data, referring donor to duplicateCorrelation; 
    ## using an adaptive span (number of genes, based on the number of genes in the dge) for smoothing the mean-variance trend
    v.swt <- voomLmFit(dge,design = des, block = as.factor(dge$samples[[batch]]), adaptive.span = T, sample.weights = T)
    
    v.swt.fit.e <- eBayes(v.swt)
    
    ## run top table over contrasts
    v.swt.e.tt <- topTable(v.swt.fit.e,coef = "APOE_carrier_E4", number=Inf, adjust.method = "BH") |>
                                 mutate(data_type = opt$datatype, 
                                        cluster = clus,
                                        cell_type_anno = c_cell_type,
                                        SpX = c_SpX,
                                        contrast = "carrier", 
                                        .before = 1) |>
                                 arrange(adj.P.Val)
    
    # names(v.swt.e.tt) <- colnames(cont)
    
    message("Done - Save data")
    saveRDS(v.swt.e.tt, file = here(data_dir, sprintf("voomLmFit_%s_%s.rds", opt$datatype, clus)))
    return(tibble(cluster = clus ,pval01 = sum(v.swt.e.tt$P.Value < 0.1), FDR05 = sum(v.swt.e.tt$adj.P.Val < 0.05)))
}, otherwise = NA)
)

write.csv(lmf_summary, file = here(data_dir, sprintf("vlmf_FDR05_summary-%s.csv", opt$datatype)), row.names = FALSE)

# slurmjobs::job_single('06_Clusterwise_voomLmFit_Xenium__cell_type_anno', create_shell = TRUE, memory = '25G', command = "Rscript 06_Clusterwise_voomLmFit_Xenium.R --datatype Xenium__cell_type_anno")
# slurmjobs::job_single('06_Clusterwise_voomLmFit_Xenium_cell_type_anno_SpX', create_shell = TRUE, memory = '25G', command = "Rscript 06_Clusterwise_voomLmFit_Xenium.R --datatype Xenium_cell_type_anno_SpX")
# slurmjobs::job_single('06_Clusterwise_voomLmFit_Xenium_Oligo.3_Astro', create_shell = TRUE, memory = '25G', command = "Rscript 06_Clusterwise_voomLmFit_Xenium.R --datatype Xenium_Oligo.3_Astro")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()