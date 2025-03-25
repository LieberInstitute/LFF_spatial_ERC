## Louise Huuki-Myers, March 2025
## Add SpD annotaions to spe, explore and annotate SpDs

library("spatialLIBD")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")
library("cowplot")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

plot_dir <- here("plots", "05_spe_correct_cluster", "19_SpD_update_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "19_SpD_update_spe")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))
spe

table(spe$BayesSpace_SVGm_k09, spe$sample_id)
table(spe$BayesSpace_SVGm_k09_anno)

## clean up colData
colnames(colData(spe))

colData(spe) <- colData(spe)[,!grepl("Markers", colnames(colData(spe)))]
colData(spe)$BayesSpace_SVGm_k09_anno <- NULL

#### load SpD annotations ####
spd_anno <- readxl::read_excel(here("processed-data","05_spe_correct_cluster", "10_spatial_registration_DLPFC", "ERC_SpD_spatial_registration_Annotations.xlsx"))

spd_anno$cluster %in% spe$BayesSpace_SVGm_k09

anno_table <- spd_anno |>
    mutate(Annotation = fct_reorder(Annotation, order)) |>
    mutate(SpD = fct_reorder(paste0(Annotation, "~", cluster), order)) |>
    select(cluster, Annotation, SpD) |>
    column_to_rownames("cluster")

all(rownames(anno_table) %in% spe$BayesSpace_SVGm_k09)

anno_table <- anno_table[as.character(spe$BayesSpace_SVGm_k09),]
dim(anno_table)

table(anno_table$SpD)

## add to colData
colData(spe) <- cbind(colData(spe), anno_table) 

table(spe$SpD, spe$BayesSpace_SVGm_k09)

#### Define colors for SpD ####

SpD_colors <- c("Vasc~Sp09D08" = "#E05AD2",
                "L1~Sp09D05" = "#021380",
                "L2.3~Sp09D01" = "#FEAF16",
                "L3~Sp09D02" = "#00BCF9",
                "L4.inhib~Sp09D09" = "#C82100",
                "L5~Sp09D03" = "#16FF32",
                "L6~Sp09D04" = "#116A52",
                "WM.uf~Sp09D07" = "#E4E1E3",
                "WM~Sp09D06" = "#500802")


color_test <- vis_clus(
    spe = spe,
    sampleid = "Br5517",
    clustervar = "SpD",
    colors = SpD_colors,
    point_size = 1.3
)
ggsave(color_test, filename = here(plot_dir, "SpD_color_test_Br5517.png"))


vis_clus_plots <- vis_grid_clus(
    spe = spe,
    clustervar = "SpD",
    colors = SpD_colors,
    sort_clust = FALSE,
    return_plots = TRUE,
    point_size = 1,
    pdf_file = here(plot_dir, "SpD_color_test.pdf")
)

apoe_anc <- as.data.frame(colData(spe)) |> 
    select(sample_id, APOE, Ancestry) |> 
    unique() |>
    mutate(apoe_anc = paste0(APOE, "_", Ancestry))

all(apoe_anc$sample_id == names(vis_clus_plots))

apoe_anc_split <- splitit(apoe_anc$apoe_anc)

pdf(here(plot_dir, "ERC_SpD_vis_clus_split.pdf"), height = 24, width = 36)
print(cowplot::plot_grid(plotlist = split_plots))
dev.off()

# pdf(here(plot_dir, "ERC_SpD_vis_clus_split.pdf"), height = 12, width = 18)
# walk2(apoe_anc_split, names(apoe_anc_split), function(index,name){
#     
#     split_plots <- vis_clus_plots[index]
#     message(name, "; ", length(split_plots))
#     plots <- cowplot::plot_grid(plotlist = split_plots)
#     
#     # title
#     title <- ggdraw() + 
#         cowplot::draw_label(
#             name,
#             fontface = 'bold',
#             x = 0,
#             hjust = 0
#         ) +
#         theme(
#             # add margin on the left of the drawing canvas,
#             # so title is aligned with left edge of first plot
#             plot.margin = margin(0, 0, 0, 7)
#         ) 
#         
#    print(cowplot::plot_grid(
#         title, plots,
#         ncol = 1,
#         # rel_heights values control vertical title margins
#         rel_heights = c(0.1, 1)
#     ))
# })
# dev.off()
# 
# cowplot_title <- function(my_title, plots){
#     title <- ggdraw() + 
#         cowplot::draw_label(
#             my_title,
#             fontface = 'bold',
#             x = 0,
#             hjust = 0
#         ) +
#         theme(
#             # add margin on the left of the drawing canvas,
#             # so title is aligned with left edge of first plot
#             plot.margin = margin(0, 0, 0, 7)
#         ) 
#     
#     print(cowplot::plot_grid(
#         title, plots,
#         ncol = 1,
#         # rel_heights values control vertical title margins
#         rel_heights = c(0.1, 1)
#     ))
# }

## plot by APOE
apoe_split <- splitit(apoe_anc$APOE)
map2(apoe_split, names(apoe_split), 
      ~vis_grid_clus(
          spe = spe[,spe$sample_id %in% apoe_anc$sample_id[.x]],
          clustervar = "SpD",
          colors = SpD_colors,
          sort_clust = FALSE,
          return_plots = FALSE,
          point_size =  2+1/length(apoe_split),
          pdf_file = here(plot_dir, sprintf("ERC_SpD_%s.pdf", gsub("/","",.y)))
      )
      )

## save all color output
save(SpD_colors, file = here("processed-data", "05_spe_correct_cluster", "SpD_colors.Rdata"))

#### Plot reduced dims ####
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(spe,
                          prefix = "ERC_spe",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "SpD",
                          color_pal = SpD_colors))



#### summary by cluster ####
message(Sys.time(), " - Summarize cluster")

pd <- as.data.frame(colData(spe))

cluster_info <- pd |>
    group_by(BayesSpace_SVGm_k09, SpD, Annotation) |>
    summarize(n = n(),
              prop = n/ncol(spe),
              median_sum_umi = median(sum_umi),
              median_sum_gene = median(sum_gene),
              median_chrM_ratio = median(expr_chrM_ratio),
              median_Nmask = median(Nmask_dark_blue)
    )

write.csv(cluster_info, file = here(data_dir, "ERC_spe_cluster_info.csv"), row.names = FALSE)

cluster_metrics_long <- pd |> 
    select(SpD, sum_gene, sum_umi, expr_chrM_ratio)  |>
    pivot_longer(!c(SpD), names_to = "metric") |>
    group_by(SpD, metric) |>
    summarize(median = median(value))

qc_violin_plot_all <- pd |> 
    select(SpD,
           sum_gene, 
           sum_umi, 
           expr_chrM_ratio)   |>
    pivot_longer(!c(SpD), names_to = "metric") |>
    ggplot() +
    geom_violin(aes(x = SpD, y = value, fill = SpD), 
                scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
    scale_fill_manual(values = SpD_colors, guide = "none") +
    theme_bw() +
    facet_grid(metric~., scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Quality Metric Value", x = "SpD")

ggsave(qc_violin_plot_all, filename = here(plot_dir, "ERC_SpD_QCmetricViolin_ALL.png"), width = 12, height = 12)


####  SpD Proportions ####

SpD_GM_proportions <- pd |>
    group_by(sample_id, APOE, Sex, Ancestry, GM = grepl("L", SpD)) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n))

sample_GM_order <- SpD_GM_proportions |> filter(GM) |> arrange(prop) |> pull(sample_id)

SpD_proportions <- pd |>
    group_by(sample_id, APOE, Sex, Ancestry, SpD) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = sample_GM_order))

SpD_proportion_bar <- SpD_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = SpD)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = SpD_colors, guide = "none") +
    labs(y = "SpD Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(SpD_proportion_bar, filename = here(plot_dir, "ERC_sn_barplot_SpD_prop.png"), width = 10)

SpD_proportion_bar_APOE <- SpD_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = SpD)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = SpD_colors) +
    facet_grid(.~APOE, scales = "free_x", space = "free") +
    labs(y = "SpD Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(SpD_proportion_bar_APOE, filename = here(plot_dir, "ERC_sn_barplot_SpD_prop_APOE.png"), width = 10)

## evaluate proportion sample by cell 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
# we need sample_colors

sample_proportions <- pd |>
    group_by(sample_id, SpD) |>
    summarize(n = n()) |> 
    group_by(SpD) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = names(sample_colors))) 

## check for SpD from few samples 
sample_proportions |> arrange(-prop)

sample_qc <- sample_proportions |> 
    group_by(SpD) |>
    summarize(n_samples = length(unique(sample_id)),
              max_prop = max(prop)) |>
    mutate(pass_sampleQC = max_prop < 0.5 & n_samples > 5)

## lowest contibution is 19 samples 
sample_qc |> arrange(n_samples)
sample_qc |> arrange(-max_prop)

## plot sample by SpD proportions bar
sample_proportion_bar <- sample_proportions |>
    ggplot(aes(x = SpD, y = prop, fill = sample_id)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > 0.3, sprintf("%s - %g", sample_id, round(prop, 3)), "")),
              position = position_stack(vjust = .5),
              size = 2, angle = 90) +
    theme_bw() +
    scale_fill_manual(values = sample_colors) +
    labs(y = "SpD Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(sample_proportion_bar, filename = here(plot_dir, "ERC_SpD_barplot_sample_prop.png"), width = 10)



#### Save data ####
message(Sys.time(), " - Saving HDF5 SPE")
saveHDF5SummarizedExperiment(
    spe,
    dir = here("processed-data", "spe_objects", "spe_ERC_annotated"),
    replace = TRUE
)

# slurmjobs::job_single('19_SpD_update_spe', create_shell = TRUE, memory = '25G', command = "Rscript 19_SpD_update_spe.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
