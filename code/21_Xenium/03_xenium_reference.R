## Louise Huuki-Myers, Feb 2026
## Create Reference file for Xenium probe design
## referencing https://github.com/LieberInstitute/xenium_NAC/blob/291f7c0174e19b902b1d563f58ed11d01849f8cb/code/00_Xenium_Panel_Design/GenerateXeniumPanelMexFiles.R

#### Set up ####
library("here")
library("sessioninfo")
library("spatialLIBD")
library("DropletUtils")
library("readxl")

data_dir <- here("processed-data", "21_Xenium", "03_xenium_reference")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


#### load sce data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

colnames(sce) <- sce$Barcode

#the recommendation via 10x is ~50k cells.
50000/ncol(sce) #[1] 50k/122004 = ~40%

cell_types <- sce$cell_type_anno

# Split cells by type
split_cells <- split(seq_along(cell_types), cell_types)

# Sample 50% from each group
set.seed(223) # For reproducibility
sampled_cells <- unlist(lapply(split_cells, function(cells) {
    sample(cells, size = length(cells) * 0.30)
}))

length(sampled_cells) #36k

sce <- sce[, sampled_cells]

table(sce$cell_type_anno)
# Astro.1          Astro.2          Astro.3          Astro.4          Astro.5            Macro          Micro.1 
# 3292             1582              926              339              243              129             1268 
# Micro.2          Micro.3          Micro.4          Micro.5            OPC.1            OPC.2            OPC.3 
# 1124              562              518              177             1444             1202             1179 
# OPC.4            OPC.5          Oligo.1          Oligo.2          Oligo.3          Oligo.4          Oligo.5 
# 588              137             3755             3603             2576             2073             1276 
# Vasc.Endo          Vasc.PC        Vasc.VLMC         Excit.L2     Excit.L2_5.1     Excit.L2_5.2       Excit.L5.1 
# 175               44              274              662             2601              182              670 
# Excit.L5.2    Excit.L5_6_NP      Excit.L6_CT        Excit.L6b Inhib.Chandelier Inhib.Lamp5_Lhx6       Inhib.Pax6 
# 72               45              187              269               90              370               37 
# Inhib.Pvalb        Inhib.Sst        Inhib.Vip 
# 408              244             2260

#are there any duplicated barcodes? 
any(duplicated(colnames(sce)))

## convert to sparse Matrix
class(counts(sce)) #DelayedMatrix

message(Sys.time(), " - convert to Sparse Matrix")
counts(sce) <- as(counts(sce), "sparseMatrix")


# Use Droplet.Utils write10xCounts() to make the mex files. 
message(Sys.time(), " - Write Mex file")
write10xCounts(path = here(data_dir, "Panel_Design_Files"),
               x = counts(sce),
               gene.id = rowData(sce)$gene_id,
               gene.symbol = rowData(sce)$gene_name,
               barcodes = colnames(sce),
               type = "sparse",
               version = "3")

list.files(here(data_dir, "Panel_Design_Files"))
# [1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"

## From https://github.com/LieberInstitute/xenium_NAC/blob/291f7c0174e19b902b1d563f58ed11d01849f8cb/code/00_Xenium_Panel_Design/GenerateXeniumPanelMexFiles.R
## &  https://github.com/LieberInstitute/spatialAmygdala/blob/devel/code/09_xenium_panel/format_for%20_10x_submission.R
# ====== Adding spatial domain annotations ======
# bundleOutputs is a function provided by 10x. Note that their code had a typo where
# the "data" argument was not actually called and instead directly used "seurat_obj" corrected below.

# Define function
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
    
    if (require("data.table", quietly = TRUE)) {
        data.table::fwrite(
            data.table::data.table(
                barcode = barcodes,
                annotation = unlist(data[[cell_type]])
            )[subset, ],
            file.path(out_dir, "annotations.csv")
        )
    } else {
        write.table(
            data.frame(
                barcode = barcodes,
                annotation = unlist(data[[cell_type]])
            )[subset, ],
            file.path(out_dir, "annotations.csv"),
            sep = ",", row.names = FALSE
        )
    }
    
    bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))
    
    utils::zip(
        bundle,
        list.files(out_dir, full.names = TRUE),
        zip = "zip"
    )
    
    if (file.info(bundle)$size / 1e6 > 500) {
        warning("The output file is more than 500 MB and will need to be subset further.")
    }
}

bundleOutputs(out_dir = here(data_dir, "Panel_Design_Files"), 
              data = sce, 
              barcodes = colnames(sce),
              cell_type = "cell_type_anno")

list.files(here(data_dir, "Panel_Design_Files"))
# [1] "annotations.csv"        "barcodes.tsv.gz"        "features.tsv.gz"        "matrix.mtx.gz"         
# [5] "Panel_Design_Files.zip"

# read annotations.csv to check
test <- read.csv(here(data_dir, "Panel_Design_Files","annotations.csv"))
head(test)

# barcode annotation
# 1 TCATGTTCACAATTCG-1_Br5161    Astro.1
# 2 CCGGGTAAGCCGCTTG-1_Br5529    Astro.1
# 3 GTGATGTGTGCGGATA-1_Br5415    Astro.1
# 4 GACTATGGTCACCGCA-1_Br6476    Astro.1
# 5 TGTGAGTAGTTTGTCG-1_Br5212    Astro.1
# 6 GATCATGGTATCACGT-1_Br1556    Astro.1

table(test$annotation)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()