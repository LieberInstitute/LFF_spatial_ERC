library(basilisk)
library(SingleCellExperiment)
library(SpatialExperiment)
library(HDF5Array)
library(zellkonverter)
library(sessioninfo)
library(here)

task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
stopifnot(task_id %in% c(1, 2))

if (task_id == 1) {
    sce_dir = here('processed-data', 'sce_objects', 'sce_ERC_subcluster')
    out_path = here('processed-data', '17_cell_cell_comm', 'sce.h5ad')
} else {
    sce_dir = here('processed-data', 'spe_objects', 'spe_ERC_annotated')
    out_path = here('processed-data', '17_cell_cell_comm', 'spe.h5ad')
}

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

#   Discard parts of objects we don't need, ensure proper data types, and
#   otherwise prepare objects for conversion
if (task_id == 1) {
    assays(sce) = list(
        counts = as(assays(sce)$counts, "dgCMatrix"),
        logcounts = as(assays(sce)$logcounts, "dgCMatrix")
    )
    reducedDims(sce) = list(UMAP = reducedDims(sce)$UMAP)
    metadata(sce) = list(cell_type_colors = metadata(sce)$cell_type_colors)
} else {
    reducedDims(sce) = list(UMAP = reducedDims(sce)$UMAP)
    
    #   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
    #   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
    #   resulting AnnData, which corresponds to reducedDims(spe)$spatial in R
    reducedDims(sce)$spatial = spatialCoords(sce)
}

write_anndata(sce, out_path)

session_info()
