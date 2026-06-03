## Louise Huuki-Myers June 2026
## Run voomLmFit on Xenium Oligo.3 by SPX

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

opt <- list()
opt$datatype <- "Xenium"

message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

spe <- spe[,spe$spot_class == "singlet"]

spe

message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, basename(pb_fn)))

#### Set up dirs ####
data_dir <- here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", sprintf("vlmf_%s", opt$datatype))
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
sce_pb <- readRDS(pb_fn)

dim(sce_pb)
table(sce_pb$registration_variable)

sce_pb$APOE_carrier_syn <- gsub("\\+", "", sce_pb$APOE_carrier)
table(sce_pb$APOE_carrier_syn)

clusters <- levels(sce_pb$registration_variable)
names(clusters) <- clusters

message(Sys.time(), " - Loop voomlmFit by cluster")

lmf_summary <- map_dfr(clusters, possibly(function(clus){
    
    dge <- sce_pb[,sce_pb$registration_variable ==clus]

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
                                        contrast = "carrier", 
                                        .before = 1) |>
                                 arrange(adj.P.Val)
    
    # names(v.swt.e.tt) <- colnames(cont)
    
    message("Done - Save data")
    saveRDS(v.swt.e.tt, file = here(data_dir, sprintf("voomLmFit_%s_%s.rds", opt$datatype, clus)))
    return(tibble(clsuter = clus ,pval05 = sum(v.swt.e.tt$P.Value < 0.05), FDR05 = sum(v.swt.e.tt$adj.P.Val < 0.05)))
}, otherwise = NA)
)

write.csv(lmf_summary, file = here(data_dir, sprintf("vlmf_FDR05_summary-%s.csv", opt$datatype)), row.names = FALSE)

# slurmjobs::job_single('01_Clusterwise_voomLmFit_Xenium', create_shell = TRUE, memory = '25G', command = "Rscript 01_Clusterwise_voomLmFit.R --datatype Xenium")

# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "01_Clusterwise_voomLmFit",
#                     create_script = FALSE)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()