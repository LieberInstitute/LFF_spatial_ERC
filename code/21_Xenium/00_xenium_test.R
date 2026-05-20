library("here")
library("spatialLIBD")
library("SpatialExperiment")

data_dir <- here("processed-data", "21_Xenium", "00_xenium_test")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# library("XeniumIO")

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

lobstr::obj_size(spe) #40.24 MB

qs2::qs_save(spe, file = here(data_dir, "spe_test_Br1039.qs2"))


#### subset test ####
spe <- qs2::qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

spe <- spe[,spe$BrNum %in% c("Br1039", "Br1556")]
lobstr::obj_size(spe) #216M

colnames(rowData(spe))[1:2] <- c('gene_id', "gene_name")
rowData(spe)$gene_search <- paste0(rowData(spe)$gene_name, "; ", rowData(spe)$gene_id)

rownames(spe) <- rowData(spe)$gene_id

## Check gene data
stopifnot(all(
    c("gene_id", "gene_name", "gene_search") %in%
        colnames(rowData(spe))
))
stopifnot(identical(
    paste0(rowData(spe)$gene_name, "; ", rowData(spe)$gene_id),
    rowData(spe)$gene_search
))

## Rownames are Ensembl ids
stopifnot(identical(rownames(spe), rowData(spe)$gene_id))

spe$key <- colnames(spe)
spe$ManualAnnotation <- spe$SpX

variables = c(
   "sum_gex",  # "sum_umi",
   "detected_gex" #"sum_gene",
)

vars_to_check <- c(
    "sample_id",
    "key",
    "ManualAnnotation",
    variables
)
if (!all(vars_to_check %in% colnames(colData(spe)))) {
    vars_missing <- vars_to_check[!vars_to_check %in% colnames(colData(spe))]
    stop(
        "Missing in colData(spe) the following variables: ",
        paste(vars_missing, collapse = ", "),
        call. = FALSE
    )
}

## A unique spot-level ID stored under spe$key
## The 'key' column is necessary for the plotly code to work.
stopifnot(length(unique(spe$key)) == ncol(spe))

## None of the values of rowData(spe)$gene_search should be re-used in
## colData(spe)
stopifnot(!all(rowData(spe)$gene_search %in% colnames(colData(spe))))

## No column named COUNT as that's used by vis_gene_p() and related
## functions.
stopifnot(!"COUNT" %in% colnames(colData(spe)))

## The counts and logcounts assays
stopifnot(length(assayNames(spe)) >= 1)

# qs2::qs_save(spe, file = here(data_dir, "spe_Xenium_test_Br1039.qs2"))
qs2::qs_save(spe, file = here(data_dir, "spe_Xenium_test.qs2"))
saveRDS(spe, file = here(data_dir, "spe_Xenium_test.rds"))



