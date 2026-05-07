library("here")
library("spatialLIBD")

data_dir <- here("processed-data", "21_Xenium", "00_xenium_test")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# library("XeniumIO")
library("here")

#### XeniumIO example ####

# zipfile <- paste0(
#     "https://mghp.osn.xsede.org/bir190004-bucket01/BiocXenDemo/",
#     "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
# )
# destfile <- XeniumIO:::.cache_url_file(zipfile)
# 

xenium_files <- list.files(here("raw-data", "xenium"), full.names = TRUE)
message("Samples: ", length(xenium_files))

list.files(here(xenium_files[[1]]), pattern = "cell_feature_matrix.tar.gz", recursive = TRUE)

## wrong format for XeniumIO
list.files(here(xenium_files[[1]], "cell_feature_matrix"))
# [1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"  

spe <- TENxXenium(resources = here(xenium_files[[1]], "cell_feature_matrix", "features.tsv.gz"),
                  xeniumOut = xenium_files[[1]],
                  sample_id = "Br1039") |>
    import(ref = "Gene Expression")

# Error in h(simpleError(msg, call)) : 
#     error in evaluating the argument 'con' in selecting a method for function 'import': The 'cell_feature_matrix.tar.gz' file

#### SpatialExperimentIO ####

library("SpatialExperimentIO")

spe <- readXeniumSXE(xenium_files[[1]])
colData(spe)
spe$sample_id <- "Br1039"

rowData(spe)

lobstr::obj_size(spe) #40.24 MB

qs2::qs_save(spe, file = here(data_dir, "spe_test_Br1039.qs2"))


#### subset test ####
spe <- qs2::qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

spe$sample_id <- gsub("01\\.", "", spe$sample_id)

spe <- spe[,spe$BrNum %in% c("Br1039", "Br1556")]
lobstr::obj_size(spe) #216M

# qs2::qs_save(spe, file = here(data_dir, "spe_Xenium_test_Br1039.qs2"))
qs2::qs_save(spe, file = here(data_dir, "spe_Xenium_test.qs2"))
saveRDS(spe, file = here(data_dir, "spe_Xenium_test.rds"))



