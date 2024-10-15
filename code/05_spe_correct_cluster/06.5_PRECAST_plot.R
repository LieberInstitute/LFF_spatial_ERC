## October 2024, Louise Huuki-Myers
## Plot spot plots for PRECAST output

library("HDF5Array")
library("spatialLIBD")
library("tidyverse")
library("Polychrome")
library("bluster")
library("pheatmap")
library("here")
library("sessioninfo")
        
#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "06_PRECAST")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "06.5_PRECAST_plot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

## PRECAST data
precast_output <- list.files(data_dir, "*VG_k*", full.names = TRUE)
names(precast_output) <- gsub(".csv", "", basename(precast_output))

precast_cluster <- map(precast_output, read.csv)

## key is in same order for all tables
all(map_lgl(precast_cluster, ~identical(.x$key, precast_cluster$PRECAST_HVG_k02$key)))

##combine
precast_tab <- as.data.frame(map_dfc(precast_cluster, ~.x$cluster))
rownames(precast_tab) <- precast_cluster$PRECAST_HVG_k02$key
# precast_tab <- bind_cols(precast_cluster$k02[,c("key"),drop = FALSE], precast_tab)
head(precast_tab)

## some spots are missing
table(spe$key %in% rownames(precast_tab))
# FALSE   TRUE 
# 162 122040

precast_tab <- precast_tab[spe$key, ]

write.csv(precast_tab, file = here(data_dir, "PRECAST_clusters.csv"))

##Add data to spe
colData(spe) <- cbind(colData(spe), precast_tab)

## which spots are missing?
table(is.na(spe$PRECAST_HVG_k02), spe$scran_discard)
table(is.na(spe$PRECAST_HVG_k02), spe$scran_low_lib_size)
#         TRUE  FALSE
# FALSE   5755 116285
# TRUE     162      0

table(is.na(spe$PRECAST_HVG_k02), spe$sample_id)
# Br1039 Br1289 Br1556 Br1691 Br1706 Br2305 Br2582 Br3974 Br5161 Br5212 Br5276 Br5367 Br5415 Br5426 Br5460 Br5517
# FALSE   3903   3648   3523   4161   3699   3069   2367   2690   4682   4674   3905   2828   3858   4824   4032   3820
# TRUE       3      0      1      0      0      0      0     32      0      0    106      0      0      0      0      3
# 
# Br5529 Br5599 Br5634 Br5712 Br5832 Br5854 Br5941 Br6085 Br6098 Br6161 Br6263 Br6321 Br6423 Br6476 Br6538
# FALSE   4556   4812   4830   3373   3546   3575   4869   4017   4569   4864   3226   3641   4315   3964   4200
# TRUE      14      0      0      0      0      0      0      0      1      0      0      0      0      0      2


#### spot plots ####
sample_ids <- sort(unique(spe$sample_id))
#   Use 'vis_grid_clus' to preserve all spots (including overlaps)

walk(colnames(precast_tab)[36:55], function(precast_name){
    
    k = parse_number(precast_name)
    
    ## create color pallet
    # cols <- Polychrome::palette36.colors(k)
    # if(k ==2) cols <- cols[seq(2)] ## fix return 3 colors bug
    # names(cols) <- sort(unique(spe[[precast_name]]))
    
    message(precast_name)
    
    plot_dir2 <- here(plot_dir, precast_name)
    if(!dir.exists(plot_dir2)) dir.create(plot_dir2, recursive = TRUE)
    
    p_list = vis_grid_clus(
        spe = spe,
        clustervar = precast_name,
        pdf_file = here(plot_dir2, paste0(precast_name, "-ALL.pdf")),
        sort_clust = FALSE,
        # colors = cols,
        spatial = FALSE,
        point_size = 1,
        sample_order = sample_ids
    )
    
    #  Vis_clus for each sample
    walk(c("Br5517"), function(samp){
        spot_plot <- vis_clus(
            spe = spe,
            sampleid = samp,
            clustervar = precast_name,
            # sort_clust = FALSE,
            # colors = cols,
            point_size = 2
        )
        ggsave(spot_plot, filename = here(plot_dir2, paste0(precast_name, "-", samp, ".pdf")))
        ggsave(spot_plot, filename = here(plot_dir2, paste0(precast_name, "-", samp, ".png")))
    })
    
})

#### how does HVG vs. SVG clustering compare? ####
library(bluster)
table(precast_tab$PRECAST_HVG_k02, precast_tab$PRECAST_SVG_k02)
table(precast_tab$PRECAST_HVG_k03, precast_tab$PRECAST_SVG_k03)

## Adjusted Rand Index
# 0.5 corresponds to “good” similarity
pairwiseRand <- map_dbl(2:28, ~pairwiseRand(precast_tab[[sprintf("PRECAST_HVG_k%02d",.x)]], precast_tab[[sprintf("PRECAST_SVG_k%02d",.x)]], mode = "index"))

compare_cluster <- tibble(k = 2:28, pairwiseRand)

plot_dir2 <- here(plot_dir, "00_explore")
if(!dir.exists(plot_dir2)) create_dir(plot_dir2)

rand_col <- ggplot(compare_cluster, aes(x = k, y = pairwiseRand)) +
    geom_col() +
    geom_hline(yintercept = 0.5)

ggsave(rand_col, file = here(plot_dir2, "pairwiseRand_col-PRECAST_HVG_vs_SVG.png"))


walk(2:28, function(k){
    message(k)
    jacc.mat <- linkClustersMatrix(precast_tab[[sprintf("PRECAST_HVG_k%02d",k)]], 
                                   precast_tab[[sprintf("PRECAST_SVG_k%02d", k)]])
    
    png(here(plot_dir2, sprintf("jacc_heatmap-PRECAST_HVG_vs_SVG_k%02d.png", k)), 
        height = 800, width = 800)
    pheatmap(jacc.mat)
    dev.off()
})




# slurmjobs::job_single('06.5_PRECAST_plot', create_shell = TRUE, memory = '10G', command = "Rscript 06.5_PRECAST_plot.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
