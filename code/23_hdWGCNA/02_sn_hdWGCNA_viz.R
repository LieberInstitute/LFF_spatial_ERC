## Louise Huuki-Myers, April 2026
## run hdWGCNA on Oligo subtypes - pick up at https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html#basic-visualization

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
library("UCell")

library("SingleCellExperiment")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "23_hdWGCNA", "02_sn_hdWGCNA_viz")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "23_hdWGCNA", "02_sn_hdWGCNA_viz")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)
# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

message(Sys.time(), " - Load data data")
seurat_obj <- readRDS(file=here("processed-data", "23_hdWGCNA", "01_sn_hdWGCNA", 'hdWGCNA_object.rds'))


###

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Oligo-M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

pdf(here(plot_dir, "KME_plots.pdf"))
print(p)
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj) |> subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])


#### Compute hub gene signature scores ####
message(Sys.time(), " - ModuleExprScor")
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

#### Module Feature Plots ####
message(Sys.time(), " - ModuleFeaturePlot")
# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)


# make a featureplot of hub scores for each module
plot_list_score <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
pdf(here(plot_dir, "ModuleFeaturePlots.pdf"))
print(wrap_plots(plot_list, ncol=6))
print(wrap_plots(plot_list_score, ncol=6))
dev.off()

#### Radar plot ####
message(Sys.time(), " - Radar Plot")
# seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$annotation), ' '))[,1]
seurat_obj$cluster <- as.character(seurat_obj$cell_type_anno)

table(seurat_obj$cluster)

pdf(here(plot_dir, "ModuleRadarPlot.pdf"))
ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster',
  barcodes = seurat_obj@meta.data |> subset(cell_type_broad == 'Oligo') |> rownames(),
  axis.label.size=4,
  grid.label.size=4
)
dev.off()


#### Network Visulization tutorial ####
# https://smorabit.github.io/hdWGCNA/articles/network_visualizations.html

message(Sys.time(), " - ModuleNetworkPlot")
ModuleNetworkPlot(
  seurat_obj,
  outdir = plot_dir
)

# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)

#### Combined hub gene network plots ####
message(Sys.time(), " - HubGeneNetworkPlot")
# hubgene network
pdf(here(plot_dir, "HubGeneNetworkPlot.pdf"))
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)

# get the list of modules:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# hubgene network
pdf(here(plot_dir, "HubGeneNetworkPlot_select.pdf"))
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = mods[1:5] # only select 5 modules
)
dev.off()

#### Applying UMAP to co-expression networks ####
message(Sys.time(), " - RunModuleUMAP")

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

message(Sys.time(), " - Plot ModuleUMAP")
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
module_umap_plot <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ggsave(module_umap_plot, filename = here(plot_dir, "module_umap_plot.png"))

pdf(here(plot_dir, "module_umap_plot_gene.png"))
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()

# slurmjobs::job_single(name = "02_sn_hdWGCNA_viz", create_shell = TRUE, memory = "50G", command = "Rscript 02_sn_hdWGCNA_viz.R")

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



