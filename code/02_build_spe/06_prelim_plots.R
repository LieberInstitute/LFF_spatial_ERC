####
## Louise Huuki-Myers Aug 2024
## create spot plots and violin plots for LFF progress report
#### 

library("SpatialExperiment")
library("spatialLIBD")
library("tidyverse")
library("DeconvoBuddies")
library("patchwork")
library("here")
library("sessioninfo")


plot_dir <- here("plots", "02_build_spe", "06_prelim_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


spe <- readRDS(here("processed-data", "02_build_spe", "spe.rds"))
dim(spe)

focus_samples <- c("Br5415", "Br5517")
layer_genes <- c("MBP", "SNAP25", "RELN", "PCP4")

## use gene names
rowData(spe)$gene_search <- rowData(spe)$gene_name
rownames(spe) <- rowData(spe)$gene_name

walk(focus_samples, function(samp){
    walk(layer_genes, function(gene){
        
        vis_gene_plot <- vis_gene(
            spe = spe,
            sampleid = samp,
            geneid = gene,
            assayname = "counts",
            point_size = 1.7
        )
        
        ggsave(vis_gene_plot, filename = here(plot_dir, paste0("vis_gene_",gene,"_",samp,".png")))
        
    })
})


colors = c(
    "#b2df8a",
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#ff7f00",
    "gold",
    "#a65628",
    "#999999",
    "black"
)

colnames(colData(spe)) <- gsub("^10x_", "", colnames(colData(spe)))

spe$kmeans_9_clusters <- paste0("k", str_padspe$kmeans_9_clusters)

test <- plot_gene_express(sce = spe[,spe$sample_id == "Br5415"],
                  genes = "MBP",
                  assay_name = "counts",
                  cat = "kmeans_9_clusters")

walk(focus_samples, function(samp){
        vis_clus_plot <- vis_clus(
            spe = spe,
            sampleid = samp,
            clustervar = "kmeans_9_clusters",
            point_size = 1.7
        )
        
        ggsave(vis_clus_plot, filename = here(plot_dir, paste0("vis_clus_k9_",samp,".png")))
})

walk(focus_samples, function(samp){
        gene_expres_plot_big <- plot_gene_express(sce = spe[,spe$sample_id == samp],
                                           genes = c("MBP", "SNAP25"),
                                           assay_name = "counts",
                                           cat = "kmeans_9_clusters", 
                                           color_pal = colors,
                                           plot_points = TRUE)        
        
        gene_expres_plot_little <- plot_gene_express(sce = spe[,spe$sample_id == samp],
                                           genes = c("RELN", "PCP4"),
                                           assay_name = "counts",
                                           cat = "kmeans_9_clusters", 
                                           color_pal = colors,
                                           plot_points = TRUE)
        
        ggsave(gene_expres_plot_big/gene_expres_plot_little, filename = here(plot_dir, paste0("gene_expres_",samp,".png")))

})


