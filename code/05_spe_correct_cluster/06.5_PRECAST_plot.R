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
data_dir <- here("processed-data", "05_spe_correct_cluster", "06_PRECAST", "00_cluster_tab")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "06.5_PRECAST_plot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the spe data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

# clear any PRECAST columns
colnames(colData(spe))
colData(spe) <- colData(spe)[,!grepl("PRECAST", colnames(colData(spe)))]

#### load PRECAST data ####
genesets <- c("HVG", "SVG", "SVGplusMarkers")
names(genesets) <- genesets

# Read in PRECAST outputs and add SpD formatting
precast_cluster <- map(genesets, function(geneset){
    ## list output files
    precast_output <- list.files(data_dir, sprintf("*%s_k*", geneset), full.names = TRUE)
    names(precast_output) <- gsub(".csv", "", basename(precast_output))
    
    ## read files
    message(Sys.time(), sprintf(" - Reading %s output, %d files", geneset, length(precast_output)))
    precast_cluster <- map2(precast_output, parse_number(names(precast_output)), 
                            ~read.csv(.x) |>
                                mutate(SpD = sprintf("Sp%02dD%02d", .y, cluster)))
    
    return(precast_cluster)
})

head(precast_cluster[[1]][[1]])
map(precast_cluster, ~names(.x[1:5]))

## Check key is in same order for all tables
map(precast_cluster, function(pcc){
    all(map_lgl(pcc, ~identical(.x$key, pcc[[1]]$key))) 
    })

## combine to one cluster table per geneset
precast_tab <- map(precast_cluster, function(cluster_data){
    ## merge in to one table
    precast_tab <- as.data.frame(map_dfc(cluster_data, ~.x$SpD))
    rownames(precast_tab) <- cluster_data[[1]]$key
    
    ## some spots are missing
    table(spe$key %in% rownames(precast_tab))
    # FALSE   TRUE 
    # 162 122040
    
    ## order and add rownames to match spe
    precast_tab <- precast_tab[spe$key, ]
    rownames(precast_tab) <- colnames(spe)
    
    return(precast_tab)
})

map(precast_tab, dim)
map(precast_tab, ~.x[1:5,1:5])

## export to CVS
walk2(precast_tab, names(precast_tab), ~write.csv(.x, file = here(data_dir, sprintf("PRECAST_%s_clusters.csv", .y))))

#### Add data to spe ####
precast_tab <- do.call("bind_cols", precast_tab)
colnames(precast_tab)

colData(spe) <- cbind(colData(spe), precast_tab)

## which spots are missing?
table(is.na(spe$PRECAST_HVG_k02), spe$scran_discard)
#         TRUE  FALSE
# FALSE   6467 115573
# TRUE     162      0

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

spe$PRECAST_FAIL <- is.na(spe$PRECAST_HVG_k02)

vis_grid_clus(
    spe = spe,
    clustervar = "PRECAST_FAIL",
    pdf_file = here(plot_dir, paste0("PRECAST_FAIL-ALL.pdf")),
    sort_clust = FALSE,
    # colors = cols,
    spatial = FALSE,
    point_size = 1,
    sample_order = sample_ids
)

# walk(colnames(precast_tab)[grep("SVGplusMarkers",colnames(precast_tab))], function(precast_name){
walk(colnames(precast_tab)[!grepl("SVGplusMarkers",colnames(precast_tab))], function(precast_name){
    
    k = parse_number(precast_name)
    
    ## create color pallet
    # cols <- Polychrome::palette36.colors(k)
    # if(k ==2) cols <- cols[seq(2)] ## fix return 3 colors bug
    # names(cols) <- sort(unique(spe[[precast_name]]))
    
    message(precast_name)
    
    plot_dir2 <- here(plot_dir, precast_name)
    if(!dir.exists(plot_dir2)) dir.create(plot_dir2, recursive = TRUE)
    
    # p_list = vis_grid_clus(
    #     spe = spe,
    #     clustervar = precast_name,
    #     pdf_file = here(plot_dir2, paste0(precast_name, "-ALL.pdf")),
    #     sort_clust = FALSE,
    #     # colors = cols,
    #     spatial = FALSE,
    #     point_size = 1.2,
    #     sample_order = sample_ids
    # )
    
    #  Vis_clus for each sample
    walk(c("Br5517", "Br5460"), function(samp){
        spot_plot <- vis_clus(
            spe = spe,
            sampleid = samp,
            clustervar = precast_name,
            # sort_clust = FALSE,
            # colors = cols,
            point_size = 2
        )
        ggsave(spot_plot, filename = here(plot_dir2, paste0(precast_name, "-", samp, ".pdf")), width = 7, height = 7)
        ggsave(spot_plot, filename = here(plot_dir2, paste0(precast_name, "-", samp, ".png")), width = 7, height = 7)
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
