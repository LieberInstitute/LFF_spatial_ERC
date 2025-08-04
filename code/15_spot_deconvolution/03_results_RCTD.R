## Louise Huuki-Myers, Aug 2025
## Compile, analyse, and plot RCTD spot deconvolution data

#### set up ####
library("ggplot2")
library("SpatialExperiment")
library("spacexr")
library("HDF5Array")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("tidyverse")
library("scatterpie")
library("patchwork")

# cell_type_col <- "cell_type_broad"
cell_type_col <- "cell_type_anno"

plot_dir <- here("plots", "15_spot_deconvolution", "03_results_RCTD", cell_type_col)
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "15_spot_deconvolution", "03_results_RCTD", cell_type_col)
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)

cell_type_broad_levels <- names(cell_type_colors$broad)
cell_type_broad_levels <- cell_type_broad_levels[cell_type_broad_levels != "Other"]

#### load results ####
message(Sys.time(), " - load data")
rctd_fn <- list.files(here("processed-data", "15_spot_deconvolution", "02_run_RCTD", cell_type_col), full.names = TRUE)
names(rctd_fn) <- gsub(sprintf("RCTD_%s-|.rds", cell_type_col), "", basename(rctd_fn))

rctd_data <- map(rctd_fn, function(fn){
    
    rctd_data <- readRDS(fn)
    
    return(list(coldata = colData(rctd_data), assays = assays(rctd_data)))
    
})

# rctd_data[[1]]$assays$weights[1:5,1:5]

##transpose to datatype
rctd_data <- list_transpose(rctd_data)
names(rctd_data)

## merge coldata
rcdt_colData <- do.call("rbind", rctd_data$coldata)

rcdt_colData$sample_id <- NULL

dim(rcdt_colData)
head(rcdt_colData)

## merge assays
rctd_assays <- rctd_data$assays

rctd_assays <- list("weights" = map(rctd_data$assays, ~.x$"weights"),
                    "weights_unconfident" = map(rctd_data$assays, ~.x$"weights_unconfident"),
                    "weights_full" = map(rctd_data$assays, ~.x$"weights_full"))

rctd_assays <- map(rctd_assays, ~do.call(cbind, .x))
map(rctd_assays, dim)


#### Load HD5F spe  ####
spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
message("spots failed RCDT: ", ncol(spe) - nrow(rcdt_colData)) #40

spe <- spe[,rownames(rcdt_colData)]


#### replace our assays with RCDT results ####
rcdt_spe <- SpatialExperiment(sample_id = spe$sample_id,
                  colData = cbind(colData(spe), rcdt_colData),
                  assays = rctd_assays,
                  spatialCoords = spatialCoords(spe),
                  imgData = imgData(spe)
                  )

assay(rcdt_spe, "cell_counts") <- ceiling(assay(rcdt_spe, "weights") * rcdt_spe$CNmask_dark_blue)
assay(rcdt_spe, "cell_counts")[1:5, 1:5]

message(Sys.time(), " - Save Data data")
saveRDS(spe, file = here(data_dir, sprintf("spe_RCTD-%s.rds", cell_type_col)))

#### Visulize cell type weights ####
message(Sys.time(), " - Visulization")

rep_sections_tb <- read.csv(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "rep_section.csv")) |>
    filter(rep_section)

# test_vis_gene <- vis_gene(rcdt_spe, 
#          # geneid = c("Oligo.1","Oligo.2", "Oligo.3","Oligo.4","Oligo.5"),
#          geneid = "Oligo",
#          assayname = "weights",
#          cont_colors = viridisLite::rocket(10, direction = -1))
# 
# ggsave(test_vis_gene, filename = here(plot_dir, "test_vis_gene.png"))

pdf(here(plot_dir, "vis_ct_ALL_rep_sections.pdf"), width = 18, height = 9)
map(rownames(rcdt_spe), function(ct){
    message("vis_ct: ", ct)
    ct_plots <- map(c("AA", "EA"), function(anc) {
        samples <- rep_sections_tb |> filter(Ancestry == anc) |> arrange(APOE)
        
        cluster_row_plots <- map(samples$sample_id, function(s) {
            vis_clus_plot <- vis_gene(
                spe = rcdt_spe,
                geneid = ct,
                assayname = "weights",
                point_size = 1,
                sampleid = s,
                cont_colors = viridisLite::rocket(10, direction = -1)
            ) 
            return(vis_clus_plot)
        })
        cluster_row <- Reduce("+", cluster_row_plots) + plot_layout(nrow = 1)
        return(cluster_row)
    })
    
    ct_grid <- Reduce("/", ct_plots)
    ggsave(ct_grid, filename = here(plot_dir, sprintf("vis_ct_%s_rep_sections.png", ct)), width = 18, height = 9)
    return(ct_grid)
})
dev.off()



#### long data ####

sample_info <- colData(spe)[,c("key", "BrNum","SpD", "APOE", "CNmask_dark_blue")] |>
    as.data.frame() 

rcdt_long <- assay(rcdt_spe, "weights") |>
    as.matrix() |>
    reshape2::melt() |>
    dplyr::rename(cell_type = Var1, key = Var2, weights = value) |>
    as_tibble() |>
    left_join(sample_info) 

rcdt_weights_summary <- rcdt_long |>
    group_by(cell_type, BrNum, SpD) |>
    summarise(max = max(weights),
              mean = mean(weights),
              median = median(weights),
              min = min(weights))

rcdt_weights_sum <- rcdt_weights_summary |> 
    group_by(BrNum, SpD) |> 
    summarise(sum = sum(mean),                                             
               n = n()) |>
    mutate(almost1 = 1 - sum < 0.01)

## which don't sum to 1?
rcdt_weights_sum |> ungroup() |> count(n)
rcdt_weights_sum |> ungroup() |> count(almost1)
rcdt_weights_sum |> arrange(sum)

#### summary plots ####

if(cell_type_col == "cell_type_broad"){
    
    rcdt_violin <- rcdt_long |>
        ggplot(aes(x = SpD, y = weights, fill = SpD)) +
        geom_violin() +
        facet_wrap(~cell_type, ncol = 1) +
        theme_bw() +
        scale_fill_manual(values = SpD_colors) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(rcdt_violin, filename = here(plot_dir, sprintf("rcdt_violin-%s.png", cell_type_col)), height = 10)
    
    
    rcdt_mean_boxplot <- rcdt_weights_summary |>
        ggplot(aes(x = SpD, y = mean, color = SpD)) +
        geom_boxplot() +
        facet_wrap(~cell_type, ncol = 1) +
        theme_bw() +
        scale_color_manual(values = SpD_colors) +
        labs(y = "mean(weight)") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(rcdt_mean_boxplot, filename = here(plot_dir, sprintf("rcdt_mean_boxplot-%s.png", cell_type_col)), height = 10)
    
    
} else if(cell_type_col == "cell_type_fine"){
    

    pdf(here(plot_dir, sprintf("rcdt_violin-%s.pdf", cell_type_col)), height = 10)
    map(cell_type_broad_levels, ~rcdt_long |>
            filter(grepl(.x, cell_type))|>
            ggplot(aes(x = SpD, y = weights, fill = SpD)) +
            geom_violin(scale = "width") +
            facet_grid(cell_type~.) +
            theme_bw() +
            labs(title = .x) +
            scale_fill_manual(values = SpD_colors) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    )
    dev.off()
    
    
    pdf(here(plot_dir, sprintf("rcdt_mean_boxplot-%s.pdf", cell_type_col)), height = 10)
    map(cell_type_broad_levels, ~rcdt_weights_summary |>
            filter(grepl(.x, cell_type))|>
            ggplot(aes(x = SpD, y = mean, color = SpD)) +
            geom_boxplot() +
            facet_grid(cell_type~.) +
            theme_bw() +
            labs(title = .x) +
            scale_color_manual(values = SpD_colors) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    )
    dev.off()
}


ggsave(rcdt_mean_boxplot, filename = here(plot_dir, "rcdt_mean_boxplot.png"), height = 10)

#### compare to sn prop ####

# load an cell_type_proportions
load(here("processed-data", "04_snRNA-seq", "28_subcluster_update_sce","cell_type_proportions.Rdata"), verbose = TRUE)
levels(cell_type_proportions$cell_type_anno)

if(cell_type_col == "cell_type_broad"){
    
    rcdt_sample_summary <- rcdt_long |>
        group_by(cell_type_broad = cell_type, sample_id = BrNum) |>
        summarise(rctd_prop = mean(weights))
    
    ## calc broad proportions 
    compare_proportions <- cell_type_proportions |>
        separate(cell_type_anno, into = c("cell_type_broad"), sep ="\\.", extra = "drop") |>
        mutate(cell_type_broad = factor(cell_type_broad, levels = cell_type_broad_levels)) |>
        group_by(sample_id, APOE, Sex, Age, Ancestry, cell_type_broad) |>
        summarise(n = sum(n), sn_prop = sum(prop)) |>
        group_by(sample_id) |>
        left_join(rcdt_sample_summary) 
    
    compare_proportions |>
        group_by(cell_type_broad) |>
        summarise(cor = cor(sn_prop, rctd_prop))
    
    # cell_type_broad     cor
    # <fct>             <dbl>
    # 1 Astro            0.648 
    # 2 Oligo            0.781 
    # 3 OPC              0.245 
    # 4 Macro            0.509 
    # 5 Micro            0.638 
    # 6 Vasc             0.0759
    # 7 Excit            0.512 
    # 8 Inhib           -0.0213
    
    compare_proportions_scatter <- compare_proportions |>
        ggplot(aes(x= sn_prop, y = rctd_prop)) +
        geom_point() +
        geom_abline() +
        facet_wrap(~cell_type_broad, nrow = 2) +
        coord_equal() +
        theme_bw()
    
    ggsave(compare_proportions_scatter, filename = here(plot_dir, "compare_proportions_scatter.png"), width = 10)
    
    
} else if(cell_type_col == "cell_type_fine"){ 
    
    rcdt_sample_summary <- rcdt_long |>
        group_by(cell_type, sample_id = BrNum) |>
        summarise(rctd_prop = mean(weights))
    
    compare_proportions <- cell_type_proportions |>
        select(sample_id, APOE, Sex, Age, Ancestry, cell_type = cell_type_anno, sn_prop = prop) |>
        left_join(rcdt_sample_summary) 
    
    prop_cor <- compare_proportions |>
        group_by(cell_type) |>
        summarise(cor = cor(rctd_prop, sn_prop))
    
    
    compare_proportions_scatter <- compare_proportions |>
        ggplot(aes(x= sn_prop, y = rctd_prop)) +
        geom_point() +
        geom_abline() +
        facet_wrap(~cell_type, nrow = 2) +
        coord_equal() +
        theme_bw()
    
    ggsave(compare_proportions_scatter, filename = here(plot_dir, sprintf("compare_proportions_scatter-%s.png", cell_type_col)), width = 10)
    
    
    pdf(here(plot_dir, sprintf("compare_proportions_scatter-%s.pdf", cell_type_col)))
    map(cell_type_broad_levels, ~ compare_proportions |>
            filter(grepl(.x, cell_type))|>
            ggplot(aes(x= sn_prop, y = rctd_prop)) +
            geom_point() +
            geom_abline() +
            facet_wrap(~cell_type, nrow = 2) +
            coord_equal() +
            labs(title = .x) +
            theme_bw()
    )
    dev.off()
    
}

# slurmjobs::job_loop(loops = list(cell_type_col = c("cell_type_broad", "cell_type_anno")),
#                     create_shell = TRUE,
#                     name = "03_results_RCTD",
#                     create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

