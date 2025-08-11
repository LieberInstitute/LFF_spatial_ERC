library(basilisk)
library(SingleCellExperiment)
library(HDF5Array)
library(zellkonverter)
library(sessioninfo)
library(here)

sce_dir = here('processed-data', 'sce_objects', 'sce_ERC_subcluster')
out_path = here('processed-data', '17_cell_cell_comm', 'sce.h5ad')

write_anndata <- function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library("zellkonverter")
                library("reticulate")

                # Convert SCE to AnnData:
                adata <- SCE2AnnData(sce)

                #  Write AnnData object to disk
                adata$write(filename = filename)

                return()
            },
            env = zellkonverterAnnDataEnv(),
            sce = sce,
            filename = out_path
        )
    )
}

sce = loadHDF5SummarizedExperiment(sce_dir)
write_anndata(sce, out_path)

session_info()
