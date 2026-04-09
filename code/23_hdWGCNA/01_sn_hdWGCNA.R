## Louise Huuki-Myers, April 2026
## run hdWGCNA on Oligo subtypes

#### Set Up ####

# single-cell analysis package
library("Seurat")
# plotting and data science packages
library("tidyverse")
library("cowplot")
library("patchwork")
# co-expression network analysis packages:
library("WGCNA")
library("hdWGCNA")

library("qs2")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "23_hdWGCNA", "01_sn_hdWGCNA")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "23_hdWGCNA", "01_sn_hdWGCNA")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

#### Read in Data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

# Convert to Seurat
message(Sys.time(), " - Convert data")
counts(sce) <- as(counts(sce), "CsparseMatrix")
logcounts(sce) <- as(logcounts(sce), "CsparseMatrix")

seurat_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
Assays(seurat_obj)

p <- DimPlot(seurat_obj, group.by='cell_type_broad', label=TRUE) +
  umap_theme() + 
  scale_color_manual(values = metadata(sce)$cell_type_colors$broad) +
  ggtitle('Huuki-Myers et al., ERC') + 
  NoLegend()

ggsave(p, filename = here(plot_dir, "Seurat_UMAP.png"))

#### Set up Seurat object for WGCNA ####
message(Sys.time(), " - Set up for WGCNA")

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "ERC_sn" # the name of the hdWGCNA experiment
)

#### Construct metacells ####
message(Sys.time(), " - Construct metacells")

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type_broad", "BrNum"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'HARMONY', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type_broad' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
message(Sys.time(), " - Normalize")
seurat_obj <- NormalizeMetacells(seurat_obj)

#### Co-expression network analysis - setDatExpr, SoftPower,  ConstructNetwork ####

message(Sys.time(), " - Specify the expression matrix")
# specify the expression matrix that we will use for network analysis - only include Oligo
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Oligo", # the name of the group of interest in the group.by column
  group.by='cell_type_broad', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = NULL, # using RNA assay
  layer = 'data' # using normalized data
)

message(Sys.time(), " - Test soft power")
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
ggsave(wrap_plots(plot_list, ncol=2), filename = here(plot_dir, "TestSoftPowers.png"))

power_table <- GetPowerTable(seurat_obj)
head(power_table)

write.csv(power_table, file = here(data_dir, "power_table.csv"))

## Construct co-expression network
message(Sys.time(), " - ConstructNetwork")
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 5,
  tom_name = 'Oligo', # name of the topoligical overlap matrix written to disk
  tom_outdir = here(data_dir, "TOM")
)

pdf(here(plot_dir, "hdWGCNA_Oligo_Dendrogram.pdf"))
PlotDendrogram(seurat_obj, main='Oligo hdWGCNA Dendrogram')
dev.off()

# load(here("code", "23_hdWGCNA", "TOM", "Oligo_TOM.rda"), verbose = TRUE)
# consTomDS

message(Sys.time(), " - Save data as backup pre-ME")
saveRDS(seurat_obj, file=here(data_dir, 'hdWGCNA_object.rds'))

#### Module Eigengenes and Connectivity ####
message(Sys.time(), " - Check TOM file")
## get TOM data
# TOM file was messed up 
GetNetworkData(seurat_obj, wgcna_name)$TOMFiles
seurat_obj@misc[['ERC_sn']][["wgcna_net"]][["TOMFiles"]] <- here(data_dir, "TOM", "Oligo_TOM.rda")

TOM <- GetTOM(seurat_obj)

message(Sys.time(), " - Scale data")
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

message(Sys.time(), " - Module Eigengenes")
# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="BrNum"
)

# harmonized module eigengenes:
message(Sys.time(), " - GetMEs")
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

#### Compute module connectivity ####
message(Sys.time(), " - ModuleConnectivity")

TOM <- GetTOM(seurat_obj)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type_anno', 
  group_name = 'Oligo'
)

# plot genes ranked by kME for each module
message(Sys.time(), " - plot genes ranked by kME")
p <- PlotKMEs(seurat_obj, ncol=5)

pdf(here(plot_dir, "genes_ranked_kME.pdf"))
print(p)
dev.off()


## Module assignment 
message(Sys.time(), " - Module assignment")
# get the module assignment table:
modules <- GetModules(seurat_obj) |> subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes - table of the top N hub genes sorted by kME
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

message(Sys.time(), " - Save data")
qs2::qs_save(sce, file = here(data_dir, "hdWGCNA_object.qs2"))
message(Sys.time(), " - Done save qs2")


# slurmjobs::job_single(name = "01_sn_hdWGCNA", create_shell = TRUE, memory = "50G", command = "Rscript 01_sn_hdWGCNA.R")

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



