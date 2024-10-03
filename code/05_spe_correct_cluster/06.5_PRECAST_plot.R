## October 2024, Louise Huuki-Myers
## Plot spot plots for PRECAST output

library("HDF5Array")
library("spatialLIBD")
library("tidyverse")
library("Polychrome")

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
precast_output <- list.files(data_dir, full.names = TRUE)
names(precast_output) <- gsub("PRECAST_HVG_(k[0-9]+).csv", "\\1", basename(precast_output))

precast_cluster <- map(precast_output, read.csv)
Reduce(identical, map(precast_cluster, ~.x$key))

## key is in same order for all tables
all(map_lgl(precast_cluster, ~identical(.x$key, precast_cluster$k02$key)))

##combine
precast_tab <- as.data.frame(map_dfc(precast_cluster, ~.x$cluster))
colnames(precast_tab) <- paste0("PRECAST_",colnames(precast_tab))
rownames(precast_tab) <- precast_cluster$k02$key
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

#### spot plots ####
sample_ids <- sort(unique(spe$sample_id))
#   Use 'vis_grid_clus' to preserve all spots (including overlaps)

walk(colnames(precast_tab), function(precast_name){
    
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
    walk(sample_ids, function(samp){
        spot_plot <- vis_clus(
            spe = spe,
            sampleid = samp,
            clustervar = precast_name,
            # sort_clust = FALSE,
            # colors = cols,
            point_size = 2
        )
        ggsave(spot_plot, filename = here(plot_dir2, paste0(precast_name, "-", samp, ".pdf")))
    })
    
})

# slurmjobs::job_single('06.5_PRECAST_plot', create_shell = TRUE, memory = '10G', command = "Rscript 06.5_PRECAST_plot.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
