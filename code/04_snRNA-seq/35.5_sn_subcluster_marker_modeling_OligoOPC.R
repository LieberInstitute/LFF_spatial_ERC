## Louise Huuki-Myers, January 2026
## Run modeling on Oligo subtypes & OPCs

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("DeconvoBuddies")
library("HDF5Array")
library("tidyverse")
library("getopt")
library("spatialLIBD")
library("scDotPlot")
library("ComplexHeatmap")

# Import command-line parameters
scec <- matrix(
    c("celltype", "c", "1", "character", "Name of celltype"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

celltype <- opt$celltype

# celltype = "OligoOPC"

data_dir <- here("processed-data", "04_snRNA-seq", sprintf("35.5_sn_subcluster_marker_modeling_%s", celltype))
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", sprintf("35.5_sn_subcluster_marker_modeling_%s", celltype))
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

## subset sce
sce <- sce[,sce$cell_type_broad %in% c("Oligo", "OPC")]

sce$cell_type_anno2 <- sce$cell_type_anno ## keep
sce$cell_type_anno <- as.character(sce$cell_type_anno)

sce$cell_type_anno <- factor(sce$cell_type_anno)
table(sce$cell_type_anno)

cell_type_colors <- metadata(sce)$cell_type_colors
load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

#  Oligo_OPC_colors
#  OPC.1     OPC.2     OPC.3     OPC.4     OPC.5   Oligo.1   Oligo.2   Oligo.3   Oligo.4   Oligo.5 
# "#D2B037" "#BDB76B" "#FFDB58" "#A2852D" "#DA9100" "#00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC"

## pick grouped of specific OPCs
if(celltype == "OligoOPC"){
    message("Use flattened OPC subtypes")
    ## flatten OPC subtypes
    sce$cell_type_anno[sce$cell_type_broad == "OPC"] <- "OPC"

    sce$cell_type_anno <- factor(sce$cell_type_anno, levels = c("OPC", "Oligo.3", "Oligo.4", "Oligo.5", "Oligo.1", "Oligo.2"))

    
    # Oligo_OPC_colors
    Oligo_OPC_colors2 <- Oligo_OPC_colors
    Oligo_OPC_colors <- c(Oligo_OPC_colors[grepl("Oligo", names(Oligo_OPC_colors))], c(OPC = "#D2B037"))

} else if(celltype == "OligoOPC2"){
    message("Use OPC subtypes")
    sce$cell_type_anno <- factor(sce$cell_type_anno, levels = c("OPC.3", "OPC.4", "OPC.1", "OPC.2", "OPC.5", "Oligo.3", "Oligo.4", "Oligo.5", "Oligo.1", "Oligo.2"))
    
}

table(sce$cell_type_anno)

#### AD risk gene dotplot ####

## read in risk genes from OpenTargets data
AD_risk <- read_csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    filter(symbol %in% rowData(sce)$gene_name) 


pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_AD_risk.pdf", celltype)))
sce |>
    scDotPlot(features = AD_risk$symbol,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = TRUE,
              groupLegends = FALSE)
dev.off()


#### Run Mean Ratio ####
marker_stats_fn <- here(data_dir, sprintf("marker_stats_MeanRatio_%s.Rdata", celltype))

if(file.exists(marker_stats_fn)){
    message(Sys.time(), " - Load marker stats data")
    load(marker_stats_fn, verbose = TRUE)
} else {
    message(Sys.time(), " - Run MeanRatio on ", celltype)
    marker_stats_MeanRatio <- get_mean_ratio(
        sce = sce,
        assay_name = "logcounts",
        cellType_col = cell_type_anno,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )
    
    save(marker_stats_MeanRatio, file = marker_stats_fn)
}


top_MeanRatio_genes <- marker_stats_MeanRatio |> 
    arrange(cellType.target) |>
    filter(MeanRatio > 1, MeanRatio.rank <= 10)

write.csv(top_MeanRatio_genes, here(data_dir, sprintf("subtype_MeanRatio_top10_%s.csv", celltype)), row.names = FALSE)
# top_MeanRatio_genes <- read.csv(here(data_dir, sprintf("subtype_MeanRatio_top10_%s.csv", celltype)))

#### plot top MR markers ####
message(Sys.time(), " - Plots")

## plot markers
plot_marker_express_ALL(
    sce,
    marker_stats_MeanRatio,
    pdf_fn = here(plot_dir, sprintf("sn_violin_MeanRatio_top10-%s.pdf", celltype)),
    n_genes = 10,
    rank_col = "MeanRatio.rank",
    anno_col = "MeanRatio.anno",
    gene_col = "gene",
    cellType_col = "cell_type_anno",
    color_pal = Oligo_OPC_colors,
    plot_points = FALSE
)

#### MeanRatio dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_MeanRatio_genes$cellType.target[match(rownames(sce), top_MeanRatio_genes$gene)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_MeanRatio.pdf", celltype)))
sce |>
    scDotPlot(features = top_MeanRatio_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Marker" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
sce |>
    scDotPlot(features = top_MeanRatio_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Marker" = Oligo_OPC_colors),
              clusterRows = TRUE,
              groupLegends = FALSE)

dev.off()


#### run cell type specific modeling ####
modeling_fn <- here(data_dir, sprintf("_subtype-%s.rds", celltype))
pseudobulk_fn = here(data_dir, sprintf("sce_pseudobulk_subtype-%s.rds", celltype))

remodel = FALSE

if(!remodel & file.exists(modeling_fn) & file.exists(pseudobulk_fn)){
    
    message(Sys.time(), " - Load modeling data")
    modeling_results <- readRDS(modeling_fn)
    
} else {
    
    message(Sys.time(), " - Running Modeling")
    ## make APOE syntatic
    sce$APOE <- gsub("/", "", sce$APOE)
    
     modeling_results <-registration_wrapper(
        sce = sce,
        var_registration = "cell_type_anno",
        var_sample_id = "sample_id",
        covars = c("APOE", "Sex", "Age", "Anc_Afr"),
        gene_ensembl = "gene_id",
        gene_name = "gene_name",
        min_ncells = 10,
        pseudobulk_rds_file = pseudobulk_fn
        )

    message(Sys.time(), " - Saving Data")
    saveRDS(modeling_results, file = modeling_fn)
    
}

sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_subtype-%s.rds", celltype)))

top_enrichment_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_enrichment_genes, here(data_dir, sprintf("subtype_enrichment_top10_%s.csv", celltype)), row.names = FALSE)

## all signif
top100_enrichment_genes <- sig_genes_extract(
    n = 100,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

top100_enrichment_genes <- top100_enrichment_genes |>
    # filter(fdr < 0.05) |> 
    arrange(fdr) |>
    group_by(test) |>
    mutate(rank = row_number()) |>
    arrange(test)

top100_enrichment_genes |> filter(fdr < 0.05)  |> count(test)

write.csv(top100_enrichment_genes, here(data_dir, sprintf("subtype_enrichment_top100_%s.csv", celltype)), row.names = FALSE)


#### Enrichment dot plots ####
message(Sys.time(), " - Enrichment dot plots")

rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_enrichment_genes$test[match(rownames(sce), top_enrichment_genes$gene)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_enrichment.pdf", celltype)))
sce |>
    scDotPlot(features = unique(top_enrichment_genes$gene),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Marker" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()


#### Oligo Marques markers Dotplot ####

Marques_markers <- consensus_marker_sets <- list(
    OPC = c("PTPRZ1", "PDGFRA", "SERPINE2", "CSPG4", "CSPG5", "VCAN"),
    COP = c("CD9", "NEU4", "BMP4","GPR17","VCAN"),
    NFOL = c("ARPC1B", "TMEM2", "CHN2", "MPZL1", "FRMD4A", "MOBP", "DDR1", "TSPAN2"),
    MFOL = c("CTPS1", "TMEM141", "OPALIN", "MAL", "PTGDS", "EVI2A", "EVI2B"),
    MOL = c("APOD", "SEPP1", "S100B"),
    MOL1 = c("FOSB", "DUSP1", "DNAJB1"),
    MOL2 = c("ANXA5", "KLK6", "MGST3"),
    MOL3 = c("CAR2", "CNTN2", "GAD2"),
    MOL4 = c("SERPINB1", "NEAT1"),
    MOL5 = c("CYP51A1", "DHCR24", "PDLIM2"),
    MOL6 = c("IL33", "APOE", "PTGDS")
)

Marques_markers <- map(Marques_markers, ~.x[.x %in% rownames(sce)])

plot_marker_express_List(
    sce,
    gene_list = Marques_markers,
    cellType_col = "cell_type_anno",
    pdf_fn = here(plot_dir, "sn_violin_Marques_markers.pdf"),
    color_pal = Oligo_OPC_colors
)

Marques_markers <- AnnotationDbi::unlist2(Marques_markers)

rowData(sce)$Marques_markers <- NULL
rowData(sce)$Marques_markers <- names(Marques_markers)[match(rownames(sce), Marques_markers)] 
table(rowData(sce)$Marques_markers)

Marques_palette <- c(
    OPC   = "#c27f3c",
    NFOL  = "#996EC3",
    COP   = "#33A02C",
    MFOL  = "#4cab98",
    MOL   = "#08306B",
    MOL1  = "#08519C",
    MOL2  = "#2171B5",
    MOL3  = "#4292C6",
    MOL4  = "#6BAED6",
    MOL5  = "#9ECAE1",
    MOL6  = "#C6DBEF"
)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_Marques_markers.pdf", celltype)))
sce |>
    scDotPlot(features = unique(Marques_markers),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marques_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Marques_markers" = Marques_palette),
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = unique(Marques_markers),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marques_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Marques_markers" = Marques_palette),
              clusterRows = TRUE,
              groupLegends = FALSE)
dev.off()

# VLMC genes = "TBX18", "VTN", "LUM", "COL12A"
Marques_markers2 <- list(OPC = c("PDGFRA", "CSPG4", "PTPRZ1", "PCDH15"), 
                         COP = c("VCAN",  "SOX6", "GPR17", "NEU4", "BMP4", "NKX2-2"))

Marques_markers2 <- map(Marques_markers2, ~.x[.x %in% rownames(sce)])

Marques_markers2 <- AnnotationDbi::unlist2(Marques_markers2)

rowData(sce)$Marques_markers2 <- NULL
rowData(sce)$Marques_markers2 <- names(Marques_markers2)[match(rownames(sce), Marques_markers2)] 
table(rowData(sce)$Marques_markers2)

Marques_palette2 <- c(
    OPC   = "navy",
    COP  = "pink")


pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_Marques_markers2.pdf", celltype)))
sce |>
    scDotPlot(features = unique(Marques_markers2),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marques_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Marques_markers" = Marques_palette2),
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = unique(Marques_markers2),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marques_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Marques_markers" = Marques_palette2),
              clusterRows = TRUE,
              groupLegends = FALSE)
dev.off()

#### spatial registration ####
# 
# layer_ <- list()
# 
# ## ERC
# layer_$spatialERC <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "-SpD.rds"))
# colnames(layer_$spatialERC$enrichment) <- gsub("_Sp", "~Sp", colnames(layer_$spatialERC$enrichment))
# 
# ## PsychENCODE enrichment
# layer_$snDLPFC_PEC <- list()
# layer_$snDLPFC_PEC$enrichment <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/registration_stats_LIBD.rds")
# 
# ## sestan sn enrichment
# layer_$sestan_EC <- readRDS(here("processed-data", "04_snRNA-seq", "24_external_data_check", "sestan_EC_modeling.rds"))
# 
# ## spatial HPC modeling 
# layer_$spatialHPC <- list()
# layer_$spatialHPC$enrichment <- read.csv("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)
# 
# names(layer_)
# # [1] "HumanPilot"   "spatialDLPFC" "spatialERC"   "snDLPFC_PEC"  "sestan_EC"
# 
# cor_layer <- map(layer_, function(layer_mod){
#     cor_layer <- layer_stat_cor(stats = $enrichment,
#                                  = layer_mod,
#                                 model_type = "enrichment",
#                                 top_n = 100)
#     return(cor_layer)
# })
# 
# ## Registration Annotation
# anno <- map(cor_layer, ~annotate_registered_clusters(
#     cor_stats_layer = .x,
#     confidence_threshold = 0.5,
#     cutoff_merge_ratio = 0.1
# ))
# 
# anno_summary <- map2_dfr(anno, names(anno), ~.x |> 
#                              select(cluster, layer_label) |>
#                              mutate(dataset = .y)) |>
#     pivot_wider(names_from = "dataset",
#                 values_from = "layer_label") |>
#     arrange(cluster)
# 
# ## save data
# spatial_registration <- list(cor_layer = cor_layer, anno = anno)
# save(spatial_registration, file = here(data_dir, sprintf("ERCsn_subtype_registration-%s.Rdata", celltype)))
# 
# 
# walk(names(cor_layer), function(ref){
#     message(ref)
#     pdf(here(plot_dir, sprintf("layer_stat_cor_%s_subtype-%s.pdf", ref, celltype)), 
#         width = 6 + (ncol(cor_layer[[ref]])/7))
#     print(layer_stat_cor_plot(
#         cor_stats_layer = cor_layer[[ref]],
#         # reference_colors = layer_colors[[ref]],
#         annotation = anno[[ref]]
#         # query_colors = cell_type_colors$anno
#     ))
#     dev.off()
# })

#### OligoOPC Cor Jakel2019 ####

enrichment_genes <- sig_genes_extract(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

sheets <- readxl::excel_sheets(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"))

jakel_Oligo_stats <- map_dfr(sheets[2:10], ~readxl::read_xlsx(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"), sheet = .x) |>
                                 mutate(Jakel_Oligo = paste0("Jakel_", .x)))

erc_v_jakel <- enrichment_genes |>
    inner_join(jakel_Oligo_stats, relationship = "many-to-many")

erc_v_jakel_cor <- erc_v_jakel |>
    group_by(test, Jakel_Oligo) |>
    summarise(n = n(),
              cor = cor(logFC, avg_logFC))

write_csv(erc_v_jakel_cor, file = here(data_dir, "erc_v_jakel_OligoOPC_cor.csv"))

erc_v_jakel_cor |>
    group_by(test) |> 
    slice_max(cor)

erc_v_jakel_cor |>
    group_by(Jakel_Oligo) |> 
    slice_max(cor)


(erc_v_jakel_cor_wide <- erc_v_jakel_cor |>
        select(-n) |>
        pivot_wider(names_from = "test", values_from = "cor") |>
        column_to_rownames("Jakel_Oligo") |>
        as.matrix())

pdf(here(plot_dir, "OligoOPC_Jakel_cor_logFC.pdf"), height = 4, width = 8)
Heatmap(t(erc_v_jakel_cor_wide), name = "logFC cor")
dev.off()

## jakel markers 

jakel_Oligo_markers <- readxl::read_xlsx(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"), sheet = "OL lineage_markers") |> 
    as.list()

jakel_Oligo_markers <- map(jakel_Oligo_markers, ~.x[.x %in% rownames(sce)])

names(jakel_Oligo_markers) <- gsub(" ", "", names(jakel_Oligo_markers))

jakel_palette <- c(OPCs = "#d94478",
                   `OPCs/COPs` = "#cb6b64",
                   `OPCs/COPs/imOLG` = "#d34f34",
                   `COPs` = "#be7d43",
                   `COPs/Oligo6` = "#d3a439",
                   `imOLG` = "#78722c",
                   `imOLG/COPs/Oligo6` = "#9db34c",
                   `Oligo6` = "#56b956",
                   `Oligo5/Oligo6` = "#46854c",
                   `Oligo2/Oligo6` = "#54bfab",
                   Oligo5 =  "#688ccc",
                   `Oligo1/Oligo5` = "#7268d2",
                   Oligo4 = "#9a4caa",
                   Oligo3 = "#d65ec6",
                   Oligo2 = "#d283bb",
                   Oligo1 = "#9a4669")

jakel_Oligo_markers <- jakel_Oligo_markers[names(jakel_palette)]

## Jakel Violin Plots
plot_marker_express_List(
    sce,
    gene_list = jakel_Oligo_markers,
    cellType_col = "cell_type_anno",
    pdf_fn = here(plot_dir, "sn_violin_Oligo_OPC_Jakel_markers.pdf"),
    color_pal = Oligo_OPC_colors,
    free_y = TRUE
)

jakel_marker1b <- plot_gene_express(
    sce,
    genes = c("SOX6", "BCAN", "ITPR2", "APOE", "CD74", "CDH20", "LURAP1L-AS1", "MAG", "KLK6", "OPALIN", "PLP1"),
    category = "cell_type_anno",
    color_pal = Oligo_OPC_colors,
    free_y = TRUE
)

ggsave(jakel_marker1b, filename = here(plot_dir, "sn_violin_Oligo_OPC_Jakel1b_markers.pdf"))

jakel_Oligo_markers <- AnnotationDbi::unlist2(jakel_Oligo_markers)

rowData(sce)$Jakel_markers <- NULL
rowData(sce)$Jakel_markers <- names(jakel_Oligo_markers)[match(rownames(sce), jakel_Oligo_markers)] 
table(rowData(sce)$Jakel_markers)


pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_Jakel_markers.pdf", celltype)), height = 10)
sce |>
    scDotPlot(features = unique(jakel_Oligo_markers),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Jakel_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Jakel_markers" = jakel_palette),
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = unique(jakel_Oligo_markers),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Jakel_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "Jakel_markers" = jakel_palette),
              clusterRows = TRUE,
              groupLegends = FALSE)

dev.off()

## Jakel Fig1d

#### OligoOPC Cor Green2024 ####

Green_Oligo_stats <- readxl::read_xlsx(here("external-data", "Green2024","Green2024_SuppTable2.xlsx"), sheet = "DEGs") |>
    filter(cell.type == "oligodendroglia",
!is.na(avg_log2FC)) |>
    mutate(Green_Oligo = state)

Green_Oligo_stats |> count(Green_Oligo)

erc_v_Green <- enrichment_genes |>
    inner_join(Green_Oligo_stats, relationship = "many-to-many")

erc_v_Green_cor <- erc_v_Green |>
    group_by(test, Green_Oligo) |>
    summarise(n = n(),
              cor = cor(logFC, avg_log2FC))

write_csv(erc_v_Green_cor, file = here(data_dir, "erc_v_Green2024_OligoOPC_cor.csv"))

erc_v_Green_cor |>
    group_by(test) |> 
    slice_max(cor)

erc_v_Green_cor |>
    group_by(Green_Oligo) |> 
    slice_max(cor)


(erc_v_Green_cor_wide <- erc_v_Green_cor |>
        select(-n) |>
        pivot_wider(names_from = "test", values_from = "cor") |>
        column_to_rownames("Green_Oligo") |>
        as.matrix())

pdf(here(plot_dir, "OligoOPC_Green_cor_logFC.pdf"), height = 4, width = 8)
Heatmap(t(erc_v_Green_cor_wide), name = "logFC cor")
dev.off()


####  GO analysis ####
library("org.Hs.eg.db")
library("clusterProfiler")

sig_enrichment_genes <- sig_genes_extract(
    n = nrow(modeling_results$enrichment),
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
) |> filter(fdr < 0.001)

sig_enrichment_genes |> count(test)


## ENTREZID look up
entrez_search <- bitr(modeling_results$enrichment$ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

DE_entrez <- sig_enrichment_genes |> 
    left_join(entrez_search, by = c("ensembl" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) 

DE_entrez |> count(test)

## Run GO 
universe <- unique(entrez_search$ENTREZID)
length(universe)

ont <- c("CC","BP","MF")
names(ont) <- ont

go_result <- map(ont, ~compareCluster(ENTREZID ~ test,
                                      data = DE_entrez, 
                                      OrgDb = org.Hs.eg.db,
                                      fun = enrichGO,
                                      universe = universe,
                                      ont = .x, ##ALL,CC,BP,MF
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      # qvalueCutoff = 0.05,
                                      readable = TRUE))

compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y))
compare_clus |> count(Cluster)

write.csv(compare_clus, file = here(data_dir, sprintf("GO_results_enrichment_%s.csv", celltype)), row.names = FALSE)

## GO dotplot

pdf(file = here(plot_dir, sprintf("GO_dotplot_enrichment_%s.pdf", celltype)), width = 10, height = 10)
walk2(go_result, names(go_result), 
      ~print(
          dotplot(.x, 
                  x = "Cluster", 
                  showCategory = 3, 
                  label_format = 60)  +
              ggtitle(paste("GO Enrichment:", .y)) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
      )
)
dev.off()

#### OligoOPC Cor Green2024 ####

# disease associated 
#### Cell type checks ####

## from https://www.biocompare.com/Editorial-Articles/590587-A-Guide-to-Oligodendrocyte-Markers/
lit_markers <- list(OPC = c("PDGFRA", "CSPG4", "MAG", "CNP", "A2B5"),
                    Oligo = c("PLP1", "ZFP191", "ZFP488", "ZFP536", "SOX17", "NKX6-2", "SMARCA4", "CD82", "TFR", "MAL"),
                    premyelin_Oligo = c("SOX10", "OLIG1", "OLIG2", "NKX2-2", "CD9"),
                    myelinating_Oligo = c("BMP4", "ENPP4", "ASAP", "TMEM10", "MOG"),
                    disease_associated = c("SERPINA3", "C4B", "TNFRSF1A", "IL1B", "IL33", "HMOX1", "TNF", "ERK", "ERK2"), #https://doi.org/10.1038/s41593-025-01873-x
                    AD_risk = c("APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1")
)


lit_markers <- map(lit_markers, ~.x[.x %in% rownames(sce)])

plot_marker_express_List(
    sce,
    gene_list = lit_markers,
    cellType_col = "cell_type_anno",
    pdf_fn = here(plot_dir, "sn_violin_lit_markers_OligoOPC.pdf"),
    color_pal = Oligo_OPC_colors
)

## Bernie genes of interest
oligo_bernie_genes <- plot_gene_express(
    sce,
    genes = c("PSEN1", "NCSTN"),
    category = "cell_type_anno",
    color_pal = Oligo_OPC_colors
)

ggsave(oligo_bernie_genes, filename = here(plot_dir, "sn_violin_PSEN1_NCSTN_OligoOPC_check.png"))


## oligo marker dot plot
lit_markers <- AnnotationDbi::unlist2(lit_markers)

rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- names(lit_markers)[match(rownames(sce), lit_markers)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_lit.pdf", celltype)))
sce |>
    scDotPlot(features = lit_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()


#### NE Genes ####
NE_receptor_genes <- c("ADRA1A","ADRA1B","ADRA1D","ADRA2A","ADRA2B","ADRA2C","ADRB1","ADRB2","ADRB3")
all(NE_receptor_genes %in% rownames(sce))

oligo_NE_receptor_expression <- plot_gene_express(
    sce,
    genes = NE_receptor_genes,
    category = "cell_type_anno",
    color_pal = Oligo_OPC_colors
)

ggsave(oligo_NE_receptor_expression, filename = here(plot_dir, "sn_violin_NE_receptor_gene_OligoOPC_check.png"))


## By carrier status 
sce$cell_type_carrier <- paste(sce$cell_type_anno2, sce$APOE_carrier)
table(sce$cell_type_carrier)

Oligo_carrier_colors <- Oligo_OPC_colors2[rep(names(Oligo_OPC_colors2), each = 2)]
names(Oligo_carrier_colors) <- paste(names(Oligo_carrier_colors), rep(levels(sce$APOE_carrier), length(Oligo_OPC_colors2)))

oligo_NE_receptor_expression_carrier <- plot_gene_express(
    sce[,grepl("Oligo|OPC.5", sce$cell_type_anno2)],
    genes = c("ADRA1A","ADRA1B","ADRB1"),
    category = "cell_type_carrier",
    color_pal = Oligo_carrier_colors,
    ncol = 1
) 

ggsave(oligo_NE_receptor_expression_carrier, filename = here(plot_dir, "sn_violin_NE_receptor_gene_OligoOPC_carrier.png"))


pdf(here(plot_dir, "sn_dotplot_NE_receptor_gene_OligoOPC.pdf"))
sce |>
    scDotPlot(features = NE_receptor_genes,
              group = cell_type_anno,
              groupAnno = cell_type_anno,
              scale = FALSE,
              annoColors = list(cell_type_anno = Oligo_OPC_colors),
              clusterRows = TRUE,
              groupLegends = FALSE)
sce |>
    scDotPlot(features = NE_receptor_genes,
              group = cell_type_anno,
              groupAnno = cell_type_anno,
              scale = TRUE,
              annoColors = list(cell_type_anno = Oligo_OPC_colors),
              clusterRows = TRUE,
              groupLegends = FALSE)
dev.off()


# slurmjobs::job_single('35.5_sn_subcluster_marker_modeling_OligoOPC', create_shell = TRUE, memory = '25G', command = "Rscript 35.5_sn_subcluster_marker_modeling_OligoOPC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

