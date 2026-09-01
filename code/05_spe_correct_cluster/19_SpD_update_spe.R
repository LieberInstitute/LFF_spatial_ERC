## Louise Huuki-Myers, March 2025
## Add SpD annotaions to spe, explore and annotate SpDs

library("spatialLIBD")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")
library("DeconvoBuddies")

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

table(spe$BayesSpace_SVGm_k11, spe$sample_id)

## clean up colData
colnames(colData(spe))

#### Update Ancestry data ####
## updated FLARE data from Feb 2025
load(here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.Rdata"), verbose = TRUE)

samples_ancestry <- samples_ancestry |> column_to_rownames("BrNum")

spe$Anc_Afr <- samples_ancestry[spe$sample_id, "YRI"]
spe$Anc_Eur <- samples_ancestry[spe$sample_id, "CEU"]

#### load SpD annotations (k11) ####
spd_anno <- readxl::read_excel(here("processed-data","05_spe_correct_cluster", "10_spatial_registration_DLPFC", "ERC_SpD_spatial_registration_anno_summary_k11.xlsx"))

spd_anno |> select(cluster, anno)

all(spd_anno$cluster %in% spe$BayesSpace_SVGm_k11)

## vSpD_anno = plain annotation (ex. "L1")
## vSpD      = "v"-prefixed simple annotation (ex. "vL1")
## vSpD_k11  = full "v"-prefixed k + domain label (ex. "vL1~Sp11D06")
## all three ordered by `order`
anno_table <- spd_anno |>
    mutate(
        vSpD_anno = fct_reorder(anno, order),
        vSpD = fct_reorder(paste0("v", anno), order),
        vSpD_k11 = fct_reorder(paste0("v", anno, "~", cluster), order)
    ) |>
    select(cluster, vSpD, vSpD_anno, vSpD_k11) |>
    column_to_rownames("cluster")

all(rownames(anno_table) %in% spe$BayesSpace_SVGm_k11)

anno_table <- anno_table[as.character(spe$BayesSpace_SVGm_k11),]
dim(anno_table)

table(anno_table$vSpD)
table(anno_table$vSpD_anno)
table(anno_table$vSpD_k11)

## add to colData
colData(spe) <- cbind(colData(spe), anno_table) 

table(spe$vSpD, spe$BayesSpace_SVGm_k11)

#### Define colors for SpD ####
load(here("processed-data", "SpD_colors.Rdata"))

# SpD_colors_V4 <- c(
#     "Vasc"      = "#E05AD2",
#     "L1"        = "#16C72B",
#     "L2"        = "#40DAF2",
#     "L3"        = "#889DF0",
#     "LD"        = "grey80",
#     "L5"        = "#0087F5",
#     "L6"        = "#021AB6",
#     "WMuf"      = "#F4A460",
#     "WMim"      = "#E8720C",
#     "WMd"       = "#581009",
#     "Inhib"     = "#C82100"
# )
# 
# # add v prefix
# names(SpD_colors_V4) <- paste0('v', names(SpD_colors_V4))

## test colors 
color_test <- vis_clus(
    spe = spe,
    sampleid = "Br5517",
    clustervar = "vSpD",
    colors = SpD_colors,
    point_size = 2,
    guide_point_size = 3
)
ggsave(color_test, filename = here(plot_dir, "vSpD_color_test_Br5517.png"))


vis_clus_plots <- vis_grid_clus(
    spe = spe,
    clustervar = "vSpD",
    colors = SpD_colors,
    sort_clust = FALSE,
    return_plots = TRUE,
    point_size = 1,
    pdf = here(plot_dir, "vSpD_color_test.pdf")
)

apoe_anc <- as.data.frame(colData(spe)) |> 
    select(sample_id, APOE, Ancestry) |> 
    unique() |>
    mutate(apoe_anc = paste0(APOE, "_", Ancestry))

all(apoe_anc$sample_id == names(vis_clus_plots))


## plot by APOE
apoe_split <- splitit(apoe_anc$APOE)
map2(apoe_split, names(apoe_split), 
      ~vis_grid_clus(
          spe = spe[,spe$sample_id %in% apoe_anc$sample_id[.x]],
          clustervar = "vSpD",
          colors = SpD_colors,
          sort_clust = FALSE,
          return_plots = FALSE,
          point_size =  2+1/length(apoe_split),
          pdf = here(plot_dir, sprintf("ERC_vSpD_%s.pdf", gsub("/","",.y)))
      )
      )

#### Plot reduced dims ####
walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"),
     ~my_plot_reduced_dim(spe,
                          prefix = "ERC_spe",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "vSpD",
                          color_pal = SpD_colors))

#### plot each sample ####
plot_dir_sample <- here("plots", "05_spe_correct_cluster", "19_SpD_update_spe", "vis_clus_sample")
if(!dir.exists(plot_dir_sample)) dir.create(plot_dir_sample, recursive = TRUE)

walk(apoe_anc$sample_id, function(s){
    vc <- vis_clus(
        spe = spe,
        sampleid = s,
        clustervar = "vSpD",
        colors = SpD_colors,
        point_size = 1.7,
        guide_point_size = 3
    )
    ggsave(vc, filename = here(plot_dir_sample, sprintf("ERC_vSpD_%s.png", s)))
})

#### plot marker genes ####
message(Sys.time(), " - Plot marker genes")

## read in marker genes from lit

lit_markers <- read_csv(here("processed-data","05_spe_correct_cluster", "00_lit_marker_genes_layer", "lit_layer_marker_summary.csv")) |>
    mutate(in_data = gene_name %in% rowData(spe)$gene_name) |>
    arrange(Layer)

## missing from our data
lit_markers |> filter(!in_data)
# gene_name Layer  n_studies studies      in_data
# <chr>     <chr>      <dbl> <chr>        <lgl>  
# 1 RELM      Layer1         1 ERC RNAScope FALSE 
# 2 KITL      Layer3         1 Ramsden et al. FALSE

lit_markers <- lit_markers |> filter(in_data)

lit_markers_list <- map(splitit(lit_markers$Layer), ~lit_markers$gene_name[.x])

plot_marker_express_List(spe, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ERC_vSpD_Layer_lit_markers.pdf"),
                         cellType_col = "vSpD",
                         gene_name_col = "gene_name",
                         color_pal = SpD_colors,
)


#### summary by cluster ####
message(Sys.time(), " - Summarize cluster")

pd <- as.data.frame(colData(spe))

cluster_info <- pd |>
    group_by(BayesSpace_SVGm_k11, vSpD, vSpD_k11) |>
    summarize(n = n(),
              prop = n/ncol(spe),
              median_sum_umi = median(sum_umi),
              median_sum_gene = median(sum_gene),
              median_chrM_ratio = median(expr_chrM_ratio),
              median_Nmask = median(Nmask_dark_blue)
    )

write.csv(cluster_info, file = here(data_dir, "ERC_spe_cluster_info.csv"), row.names = FALSE)

cluster_metrics_long <- pd |> 
    select(vSpD, vSpD_k11, sum_gene, sum_umi, expr_chrM_ratio)  |>
    pivot_longer(!c(vSpD, vSpD_k11), names_to = "metric") |>
    group_by(vSpD, vSpD_k11, metric) |>
    summarize(median = median(value))

## violin plot uses vSpD_anno (plain) to match SpD_colors
qc_violin_plot_all <- pd |> 
    select(vSpD,
           sum_gene, 
           sum_umi, 
           expr_chrM_ratio)   |>
    pivot_longer(!c(vSpD), names_to = "metric") |>
    ggplot() +
    geom_violin(aes(x = vSpD, y = value, fill = vSpD), 
                scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
    scale_fill_manual(values = SpD_colors, guide = "none") +
    theme_bw() +
    facet_grid(metric~., scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Quality Metric Value", x = "SpD")

ggsave(qc_violin_plot_all, filename = here(plot_dir, "ERC_vSpD_QCmetricViolin_ALL.png"), width = 12, height = 12)


####  SpD Proportions ####

SpD_GM_proportions <- pd |>
    group_by(sample_id, APOE, Sex, Ancestry, GM = grepl("L", vSpD)) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n))

sample_GM_order <- SpD_GM_proportions |> filter(GM) |> arrange(prop) |> pull(sample_id)

SpD_proportions <- pd |>
    group_by(sample_id, APOE, Sex, Ancestry, vSpD, vSpD_k11) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = sample_GM_order))

## save proportion data
save(SpD_proportions, file = here(data_dir, "SpD_proportions.Rdata"))
write.csv(SpD_proportions, file = here(data_dir, "SpD_proportions.csv"))

## proportion plots 
SpD_proportion_bar <- SpD_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = vSpD)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = SpD_colors, guide = "none") +
    labs(y = "SpD Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(SpD_proportion_bar, filename = here(plot_dir, "ERC_sn_barplot_vSpD_prop.png"), width = 10)

SpD_proportion_bar_APOE <- SpD_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = vSpD)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = SpD_colors) +
    facet_grid(.~APOE, scales = "free_x", space = "free") +
    labs(y = "SpD Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(SpD_proportion_bar_APOE, filename = here(plot_dir, "ERC_sn_barplot_vSpD_prop_APOE.png"), width = 10)

## evaluate proportion sample by cell 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
# we need sample_colors

sample_proportions <- pd |>
    group_by(sample_id, vSpD, vSpD_k11) |>
    summarize(n = n()) |> 
    group_by(vSpD, vSpD_k11) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = names(sample_colors))) 

## check for SpD from few samples 
sample_proportions |> arrange(-prop)

sample_qc <- sample_proportions |> 
    group_by(vSpD, vSpD_k11) |>
    summarize(n_samples = length(unique(sample_id)),
              max_prop = max(prop)) |>
    mutate(pass_sampleQC = max_prop < 0.5 & n_samples > 5)

## lowest contibution is 19 samples 
sample_qc |> arrange(n_samples)
sample_qc |> arrange(-max_prop)

## plot sample by SpD proportions bar (x-axis uses vSpD_anno, plain)
sample_proportion_bar <- sample_proportions |>
    ggplot(aes(x = vSpD, y = prop, fill = sample_id)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > 0.3, sprintf("%s - %g", sample_id, round(prop, 3)), "")),
              position = position_stack(vjust = .5),
              size = 2, angle = 90) +
    theme_bw() +
    scale_fill_manual(values = sample_colors) +
    labs(y = "SpD Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(sample_proportion_bar, filename = here(plot_dir, "ERC_vSpD_barplot_sample_prop.png"), width = 10)


#### GTF file ####
## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
all(rownames(spe) %in% gtf$gene_id)

match_genes <- match(rownames(spe), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(spe) <- gtf[match_genes]


#### Add colors to metadata ####
metadata(spe)$SpD_colors <- SpD_colors

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
sessioninfo::session_info()
