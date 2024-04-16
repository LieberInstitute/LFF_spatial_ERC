
library("spatialLIBD")
library("tidyverse")
library("gridExtra")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "02_build_spe", "02_explore_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load spe
spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))

#### Explore colData ####
colData(spe)

## no BrNum
colnames(colData(spe))
# [1] "sample_id"              "in_tissue"              "array_row"              "array_col"             
# [5] "10x_graphclust"         "10x_kmeans_10_clusters" "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters" 
# [9] "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters" 
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                    "sum_umi"               
# [17] "sum_gene"               "expr_chrM"              "expr_chrM_ratio"        "ManualAnnotation"      
# [21] "subject"                "age"                    "sex"                    "race"                  
# [25] "diagnosis"              "rin"                    "apoe"                   "Nmask_dark_blue"       
# [29] "Pmask_dark_blue"        "CNmask_dark_blue"       "overlaps_tissue"    

## 
table(spe$sample_id)
## all have 4992, this is not filererd for spots not in tissue
table(spe$in_tissue, spe$overlaps_tissue)
table(spe$sample_id, spe$overlaps_tissue)

length(unique(spe$subject))

table(duplicated(colnames(spe)))
# FALSE   TRUE 
# 4992 149760 

## compare with metadata

#### plot UMIs ####
spe <- spe[,spe$in_tissue]

samples <- unique(spe$subject)

umi_test <- spatialLIBD::vis_gene(spe, sampleid = "V13B23-363_A1", geneid = "sum_umi")
ggsave(umi_test, filename = here(plot_dir, "umi_test.png"))

umi_plots <- map(samples, ~vis_gene(spe, sampleid = .x, geneid = "sum_umi"))
## TODO multi plot pages
# umi_plot_pages <- marrangeGrob(umi_plots, nrow=3, ncol=2)
# # unable to start device X11cairo
# 
# ## arrangeGrob returns a [1] "gtable" "gTree"  "grob"   "gDesc" not compatable with ggsave
# class(umi_plot_pages)
# typeof(umi_plot_pages)
# 
# ggsave(
#     plot = umi_plot_pages$grob,
#     filename = here(plot_dir, "Sample_umi.pdf"), 
#     width = 8.5, height = 11.5
# )


pdf(here(plot_dir, "Sample_umi.pdf"))
umi_plots
dev.off()

# vis_grid_gene(spe, 
#               geneid = "sum_umi",
#               pdf_file = "Sample_umi_grid.pdf",
#               return_plots = FALSE)

k9_plots <- map(samples, ~vis_clus(spe, sampleid = .x, clustervar = "10x_kmeans_9_clusters"))

pdf(here(plot_dir, "Sample_kmeans_k9.pdf"))
k9_plots
dev.off()

