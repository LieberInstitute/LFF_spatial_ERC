
library("spatialLIBD")
library("here")

## local data
spe <- readRDS(here("processed-data", "spe_objects", "spe_raw.rds"))
# lobstr::obj_size(spe)
# 5.60 GB


qc_file = here("processed-data", "05_QC_anno", "spe_qc_anno_clean.csv")
if(file.exists(qc_file)){
    qc_anno_clean <- read.csv(qc_file)
    pd <- as.data.frame(colData(spe))
    
    pd <- pd |>
        left_join(qc_anno_clean |> select(qc_anno, key = spot_name)) |>
        mutate(qc_anno = factor(ifelse(is.na(qc_anno) & in_tissue, "None", qc_anno)),
               scran_qc_anno = factor(ifelse(as.logical(scran_discard), paste0("scran-", qc_anno), as.character(qc_anno)))
        )
    
    spe$qc_anno <- pd$qc_anno
    spe$scran_qc_anno <- pd$scran_qc_anno
}

## keep 10x UMAP
reducedDimNames(spe)
# [1] "10x_pca"  "10x_tsne" "10x_umap"

# reducedDims(spe)$`10x_tsne` <- NULL
# reducedDims(spe)$`10x_pca` <- NULL


## Drop images we don't really use in the app
imgData(spe) <- imgData(spe)[
    !imgData(spe)$image_id %in% c("hires", "detected", "aligned"),
]

## remove out_tissue spots?
# spe <- spe[, spe$in_tissue]

lobstr::obj_size(spe)
# 2.79 GB

saveRDS(spe, here("code", "02_build_spe", "02_QC_app", "spe_raw.rds"))
