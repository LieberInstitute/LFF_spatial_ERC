## Louise Huuki-Myers, Aug 2025
## create bed files for cluster marker genes, first step in LDSC pipline
## based on https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/01_aggregate.R
## & https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/02_specificity_score.R
## & on https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/03_bed.R

#### Set Up ####

library("tidyverse")
library("getopt")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type",
      "mode", "m", "1", "character", "specificity or enrichment"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype <- "Visium"
# opt$datatype <- "sn_broad"
# opt$mode <- "specificity"

print(opt)

## filenames for "cluster only" pseudobulk & enrichment data
if(opt$datatype == "Visium"){
    pb_fn <- here("processed-data", "05_spe_correct_cluster", "25_SpD_pseudobulk","spe_pseudobulk_only-SpD_syn.rds")
    enrich_fn <- here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds")
} else if(opt$datatype == "sn_broad"){
    pb_fn <- here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk","sce_ERC_subcluster_pseudobulk_only-cell_type_broad.rds")
    # enrich_fn <- 
}else if(opt$datatype == "sn_fine"){
    pb_fn <- here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk","sce_ERC_subcluster_pseudobulk_only-cell_type_anno.rds")
    # enrich_fn <- 
} else {
    stop("non-valid datatype")
}

data_dir <- here("processed-data", "18_LDSC", "01_prep_marker_bed", sprintf("bed_%s_%s", opt$mode, opt$datatype))
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### specificity ###
if(opt$mode == "specificity"){
    
    message(Sys.time(), sprintf(" - Specificity: Datatype = %s, loading '%s'", opt$datatype, basename(pb_fn)))
    
    sce_pb <- readRDS(pb_fn)
    dim(sce_pb)
    
    # colnames(sce_pb)
    
    ## log2-CPM values from pseudobulked data
    dat <- assay(sce_pb, "logcounts")
    
    is_top_expressed <- function(column, q = 0.9) {
        top_10_percent_threshold <- quantile(column, q)
        as.numeric(column >= top_10_percent_threshold)
    }
    
    res <- apply(dat, 2, is_top_expressed)
    
    res <- cbind("geneName" = rownames(dat), as.data.frame(res))
    
    head(res)
    
    colSums(res[-1])
    
}

#### make Bed file ####
message(Sys.time(), " - Make Bed Files")

## gene metadata
gene.anno.file <- "/dcs04/lieber/shared/statsgen/LDSC/base/gene_meta_hg19.txt"
gene.meta <- read.table(gene.anno.file,sep="\t",header=T)

modules <- colnames(res)[-1]

for(i in 1:length(modules)){
    idx <- which(colnames(res) == modules[i])
    genes <- res$geneName[res[,idx]==1]
    genes <- intersect(genes,gene.meta$Gene.stable.ID)
    bed <- gene.meta[is.element(gene.meta$Gene.stable.ID, genes),3:5]
    bed[,1] <- paste0("chr",bed[,1])
    bed[,2] <- ifelse(bed[,2] - 100000 > 0, bed[,2] - 100000, 0)
    bed[,3] <- bed[,3] + 100000
    idx <- bed[,1] !="chrX" & bed[,1] !="chrY" & bed[,1] !="chrMT"
    bed <- bed[idx,]
    filename <- paste0(modules[i],".bed")
    outfile <- here::here(data_dir, filename)
    write.table(bed,outfile,row.names=F,col.names=F,sep="\t",quote=F)
}

# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "01_prep_marker_bed",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()






