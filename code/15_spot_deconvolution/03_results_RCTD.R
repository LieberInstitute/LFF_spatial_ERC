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
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("cell_type_col", "c", "1", "character", "cell type column"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

cell_type_col <- opt$cell_type_col
# cell_type_col <- "cell_type_broad"
# cell_type_col <- "cell_type_anno"

message("Cell type col: ", cell_type_col)

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
spe_rctd <- SpatialExperiment(sample_id = spe$sample_id,
                  colData = cbind(colData(spe), rcdt_colData),
                  assays = rctd_assays,
                  spatialCoords = spatialCoords(spe),
                  imgData = imgData(spe)
                  )

#### Add updated cell counts ####
nuc_count_fn <- list.files(here("processed-data","16_nuclear_counts"), pattern = ".csv.gz", full.names = TRUE)

nuc_count <- map_dfr(nuc_count_fn, ~read.csv(gzfile(.x), row.names = 1))
dim(nuc_count)
head(nuc_count)

# match key
spe_rctd$key2 <- paste0(gsub("Br[0-9]+","",spe$key), spe$VNum)

table(spe_rctd$key2 %in% rownames(nuc_count))

spe_rctd$num_nuclei_within <- nuc_count[spe_rctd$key2, "num_nuclei_within"]
spe_rctd$num_nuclei_intersect <- nuc_count[spe_rctd$key2, "num_nuclei_intersect"]

pd <- colData(spe_rctd) |> as.data.frame()

pd |>
    group_by(SpD) |>
    summarise(median_nuc = median(CNmask_dark_blue),
              n_0_nuc = sum(CNmask_dark_blue == 0),
              p_0_nuc = n_0_nuc/n(),
              n_1_nuc = sum(CNmask_dark_blue == 1),
              p_1_nuc = n_1_nuc/n(),
              q59 = quantile(CNmask_dark_blue, 0.90),
              n_30_nuc = sum(CNmask_dark_blue > 30),
              p_30_nuc = n_30_nuc/n(),
              max = max(CNmask_dark_blue))

pd |>
    group_by(SpD) |>
    summarise(median_nuc = median(num_nuclei_within),
              n_0_nuc = sum(num_nuclei_within == 0),
              p_0_nuc = n_0_nuc/n(),
              n_1_nuc = sum(num_nuclei_within == 1),
              p_1_nuc = n_1_nuc/n(),
              q59 = quantile(num_nuclei_within, 0.90),
              n_30_nuc = sum(num_nuclei_within > 30),
              p_30_nuc = n_30_nuc/n(),
              max = max(num_nuclei_within))

# SpD           median_nuc n_0_nuc p_0_nuc   q59 n_30_nuc  p_30_nuc   max
# <fct>              <dbl>   <int>   <dbl> <dbl>    <int>     <dbl> <int>
# 1 Vasc~Sp09D08           3       0       0    10        9 0.00145      59
# 2 L1~Sp09D05             2       0       0     6        4 0.000260     36
# 3 L2.3~Sp09D01           4       0       0     8        0 0            29
# 4 LD~Sp09D02             2       0       0     6       10 0.000378     44
# 5 Inhib~Sp09D09          3       0       0     8        1 0.000145     31
# 6 L5~Sp09D03             4       0       0     9        1 0.0000861    39
# 7 L6~Sp09D04             3       0       0     8        0 0            27
# 8 WM.uf~Sp09D07          3       0       0     7        2 0.000161     39
# 9 WM~Sp09D06             4       0       0     8        0 0            26


## add cell count assays based on nuc counts
assay(spe_rctd, "cell_counts") <- assay(spe_rctd, "weights") * spe_rctd$num_nuclei_within
assay(spe_rctd, "cell_counts")[1:10, 1:5]


## determine major cell type
# which.max(assay(spe_rctd, "weights")[, 1])

major_cell_index <- apply(assay(spe_rctd, "weights"), 2, which.max)
spe_rctd$major_cell <- rownames(spe_rctd)[major_cell_index]

table(spe_rctd$major_cell)

message(Sys.time(), " - Save Data")
saveRDS(spe_rctd, file = here(data_dir, sprintf("spe_RCTD-%s.rds", cell_type_col)))

# spe_rctd <- readRDS(here(data_dir, sprintf("spe_RCTD-%s.rds", cell_type_col)))

#### Visualize nuc counts ####

n_nuclei_boxplot <- ggplot(pd, aes(x = SpD, y = num_nuclei_within, fill = SpD)) +
    geom_boxplot(draw_quantiles = c(.5)) +
    scale_fill_manual(values = SpD_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(n_nuclei_boxplot, filename = here(plot_dir, "ERC_Visium_SpD_boxplot_n_nuclei.png"), width = 7, height =4)


## vis_gene

vis_n_nuc_test <- vis_gene(
    spe = spe_rctd,
    geneid = "num_nuclei_within",
    assayname = "weights",
    point_size = 2,
    cont_colors = viridisLite::rocket(10, direction = -1)
)

ggsave(vis_n_nuc_test, filename = here(plot_dir, "vis_n_nuc_test.png"))

## rep sections
source(here("code", "15_spot_deconvolution", "vis_rep_sections.R"))

vis_rep_n_nuclei <- vis_rep_sections(spe_rctd, geneid = "num_nuclei_within", assayname = "weights")
ggsave(vis_rep_n_nuclei, filename = here(plot_dir, "vis_rep_n_nuclei.png"), width = 18, height = 9)

#### Visualize cell type weights ####
message(Sys.time(), " - Visualization")

## rep sections weights

pdf(here(plot_dir, "vis_ct_ALL_rep_sections.pdf"), width = 18, height = 9)
map(rownames(spe_rctd), function(ct){
    message("vis_ct: ", ct)
    
    ct_grid <- vis_rep_sections(spe_rctd, assayname = "weights", geneid = ct)
    
    ggsave(ct_grid, filename = here(plot_dir, sprintf("vis_ct_%s_rep_sections.png", ct)), width = 18, height = 9)
    return(ct_grid)
})
dev.off()


## rep sections cell counts
pdf(here(plot_dir, "vis_ct_ALL_rep_sections_counts.pdf"), width = 18, height = 9)
map(rownames(spe_rctd), function(ct){
    message("vis_ct: ", ct)
    
    ct_grid <- vis_rep_sections(spe_rctd, assayname = "cell_counts", geneid = ct, minCount = 1)
    
    ggsave(ct_grid, filename = here(plot_dir, sprintf("vis_ct_%s_rep_sections_counts.png", ct)), width = 18, height = 9)
    return(ct_grid)
})
dev.off()


#### major cell count plots ####

table(spe_rctd$sample_id, spe_rctd$major_cell)

SpD_major_cell <- table(spe_rctd$SpD, spe_rctd$major_cell)

SpD_major_cell_prop <- sweep(SpD_major_cell, 1, table(spe_rctd$SpD), FUN = '/')

library(ComplexHeatmap)

library(circlize)

col_fun = colorRamp2(c (0, max(SpD_major_cell)), c("white", "red"))
pdf(here(plot_dir, sprintf("rctd_major_%s_count.pdf", cell_type_col)))
Heatmap(SpD_major_cell, 
        name = "count",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = col_fun)
dev.off()



col_fun = colorRamp2(c (0, max(SpD_major_cell_prop)), c("white", "red"))

pdf(here(plot_dir, sprintf("rctd_major_%s_prop.pdf", cell_type_col)))
Heatmap(SpD_major_cell_prop, 
        name = "prop",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = col_fun)
dev.off()


#### long data ####

sample_info <- colData(spe_rctd)[,c("key", "BrNum","SpD", "APOE", "num_nuclei_within")] |>
    as.data.frame() 

rcdt_long <- assay(spe_rctd, "cell_counts") |>
    as.matrix() |>
    reshape2::melt() |>
    dplyr::rename(cell_type = Var1, key = Var2, cell_counts = value) |>
    as_tibble() |>
    left_join(sample_info) 

rcdt_cell_counts_summary <- rcdt_long |>
    group_by(cell_type, BrNum, SpD) |>
    summarise(max = max(cell_counts),
              mean = mean(cell_counts),
              median = median(cell_counts),
              min = min(cell_counts))

rcdt_cell_counts_sum <- rcdt_cell_counts_summary |> 
    group_by(BrNum, SpD) |> 
    summarise(sum = sum(mean),                                             
               n = n())

nonZero_summary <- rcdt_long |>
    group_by(cell_type, SpD, BrNum) |>
    summarise(non_zero = sum(cell_counts > 0),
              p_non_zero = non_zero/n(),
              one_plus = sum(cell_counts >= 1),
              p_one_plus = one_plus/n(),
    )

#### summary plots ####

if(cell_type_col == "cell_type_broad"){
    
    ## plot counts
    rcdt_violin <- rcdt_long |>
        ggplot(aes(x = SpD, y = cell_counts, fill = SpD)) +
        geom_violin(scale = "width") +
        facet_wrap(~cell_type, ncol = 1) +
        theme_bw() +
        scale_fill_manual(values = SpD_colors) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(rcdt_violin, filename = here(plot_dir, sprintf("rcdt_violin_cell_counts-%s.png", cell_type_col)), height = 10)
    
    rcdt_violin_onePlus <- rcdt_long |>
        filter(cell_counts > 1 & cell_counts <20) |>
        ggplot(aes(x = SpD, y = cell_counts, fill = SpD)) +
        geom_violin() +
        facet_wrap(~cell_type, ncol = 1) +
        theme_bw() +
        scale_fill_manual(values = SpD_colors) +
        labs(y = "cell_counts [1:20]") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(rcdt_violin_onePlus, filename = here(plot_dir, sprintf("rcdt_violin_cell_counts_onePlus-%s.png", cell_type_col)), height = 10)
    
    ## summary metrics
    rcdt_mean_boxplot <- rcdt_cell_counts_summary |>
        ggplot(aes(x = SpD, y = mean, color = SpD)) +
        geom_boxplot() +
        facet_wrap(~cell_type, ncol = 1) +
        theme_bw() +
        scale_color_manual(values = SpD_colors) +
        labs(y = "mean(weight)") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(rcdt_mean_boxplot, filename = here(plot_dir, sprintf("rcdt_mean_boxplot-%s.png", cell_type_col)), height = 10)
    
    rcdt_nonZero_boxplot <- nonZero_summary |>
        ggplot(aes(x = SpD, y = p_non_zero, color = SpD)) +
        geom_boxplot() +
        facet_wrap(~cell_type, ncol = 1) +
        theme_bw() +
        scale_color_manual(values = SpD_colors) +
        labs(y = "prop cell_count > 0") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(rcdt_nonZero_boxplot, filename = here(plot_dir, sprintf("rcdt_nonZero_boxplot-%s.png", cell_type_col)), height = 10)
     
    rcdt_onePlus_boxplot <- nonZero_summary |>
        ggplot(aes(x = SpD, y = p_one_plus, color = SpD)) +
        geom_boxplot() +
        facet_wrap(~cell_type, ncol = 1) +
        theme_bw() +
        scale_color_manual(values = SpD_colors) +
        labs(y = "prop cell_count >= 0") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(rcdt_onePlus_boxplot, filename = here(plot_dir, sprintf("rcdt_onePlus_boxplot-%s.png", cell_type_col)), height = 10)
    
} else if(cell_type_col == "cell_type_fine"){
    

    pdf(here(plot_dir, sprintf("rcdt_violin_cell_counts-%s.pdf", cell_type_col)), height = 10)
    map(cell_type_broad_levels, ~rcdt_long |>
            filter(grepl(.x, cell_type))|>
            ggplot(aes(x = SpD, y = cell_counts, fill = SpD)) +
            geom_violin(scale = "width") +
            facet_grid(cell_type~.) +
            theme_bw() +
            labs(title = .x) +
            scale_fill_manual(values = SpD_colors) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    )
    dev.off()
    
    pdf(here(plot_dir, sprintf("rcdt_violin_cell_counts_onePlus-%s.pdf", cell_type_col)), height = 10)
    map(cell_type_broad_levels, ~rcdt_long |>
            filter(grepl(.x, cell_type),
                   cell_counts > 1 & cell_counts <20) |>
            ggplot(aes(x = SpD, y = cell_counts, fill = SpD)) +
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
    map(cell_type_broad_levels, ~rcdt_cell_counts_summary |>
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
    
    pdf(here(plot_dir, sprintf("rcdt_nonZero_boxplot-%s.pdf", cell_type_col)), height = 10)
    map(cell_type_broad_levels, ~nonZero_summary |>
            filter(grepl(.x, cell_type))|>
            ggplot(aes(x = SpD, y = p_non_zero, color = SpD)) +
            geom_boxplot() +
            facet_grid(cell_type~.) +
            theme_bw() +
            labs(title = .x) +
            scale_color_manual(values = SpD_colors) +
            labs(y = "prop cell count > 0") +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    )
    dev.off()    
    
    pdf(here(plot_dir, sprintf("rcdt_onePlus_boxplot-%s.pdf", cell_type_col)), height = 10)
    map(cell_type_broad_levels, ~nonZero_summary |>
            filter(grepl(.x, cell_type))|>
            ggplot(aes(x = SpD, y = p_one_plus, color = SpD)) +
            geom_boxplot() +
            facet_grid(cell_type~.) +
            theme_bw() +
            labs(title = .x) +
            scale_color_manual(values = SpD_colors) +
            labs(y = "prop cell count >= 1") +
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

rcdt_sample_prop <- rcdt_long |>
    group_by(cell_type, sample_id = BrNum) |>
    summarise(n_ct = sum(cell_counts)) |>
    group_by(sample_id) |>
    mutate(rctd_prop = n_ct/sum(n_ct))

# rcdt_sample_prop |> summarise(sum(prop))

if(cell_type_col == "cell_type_broad"){
    
    ## calc broad proportions 
    compare_proportions <- cell_type_proportions |>
        separate(cell_type_anno, into = c("cell_type_broad"), sep ="\\.", extra = "drop") |>
        mutate(cell_type_broad = factor(cell_type_broad, levels = cell_type_broad_levels)) |>
        group_by(sample_id, APOE, Sex, Age, Ancestry, cell_type_broad) |>
        summarise(n = sum(n), sn_prop = sum(prop)) |>
        group_by(sample_id) |>
        left_join(rcdt_sample_prop |> rename(cell_type_broad = cell_type)) 
    
    compare_proportions |>
        group_by(cell_type_broad) |>
        summarise(cor = cor(sn_prop, rctd_prop))
    
    # cell_type_broad     cor
    # <fct>             <dbl>
    # 1 Astro            0.643 
    # 2 Oligo            0.775 
    # 3 OPC              0.290 
    # 4 Macro            0.508 
    # 5 Micro            0.633 
    # 6 Vasc             0.0545
    # 7 Excit            0.504 
    # 8 Inhib           -0.0472
    
    compare_proportions_scatter <- compare_proportions |>
        ggplot(aes(x= sn_prop, y = rctd_prop)) +
        geom_point() +
        geom_abline() +
        facet_wrap(~cell_type_broad, nrow = 2) +
        coord_equal() +
        theme_bw()
    
    ggsave(compare_proportions_scatter, filename = here(plot_dir, "compare_proportions_scatter.png"), width = 10)
    
    
} else if(cell_type_col == "cell_type_fine"){ 
    
    compare_proportions <- cell_type_proportions |>
        select(sample_id, APOE, Sex, Age, Ancestry, cell_type = cell_type_anno, sn_prop = prop) |>
        left_join(rcdt_sample_prop) 
    
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

