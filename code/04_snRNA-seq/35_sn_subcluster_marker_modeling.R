## Louise Huuki-Myers, January 2025
## Run DeconvoBuddies::get_mean_ratio on snRNA_seq data

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

# celltype = "Astro"

data_dir <- here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", celltype)
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", celltype)
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

## subset sce
stopifnot(celltype %in% sce$cell_type_broad)
sce <- sce[,sce$cell_type_broad == celltype]

sce$cell_type_anno <- droplevels(sce$cell_type_anno)
table(sce$cell_type_anno)

cell_type_colors <- metadata(sce)$cell_type_colors

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
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
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
        cellType_col = "cell_type_anno",
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
    color_pal = cell_type_colors$anno,
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
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "Marker" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()


#### run cell type specific modeling ####
modeling_fn <- here(data_dir, sprintf("modeling_results_subtype-%s.rds", celltype))
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
# top_enrichment_genes <- read.csv(here(data_dir, sprintf("subtype_enrichment_top10_%s.csv", celltype)), row.names = 1)

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

## Summary table

marker_summary <- top_enrichment_genes |>
    filter(top <= 5) |>
    group_by(cell_type = test) |>
    summarise(enrichment_top5 = paste(gene, collapse = ", ")) |>
    left_join(top_MeanRatio_genes |>
                  filter(MeanRatio.rank <= 5, MeanRatio > 1) |>
                  group_by(cell_type = cellType.target) |>
                  summarise(MeanRatio_top5 = paste(gene, collapse = ", "))) |>
    left_join(top_enrichment_genes |> select(cell_type = test, gene) |>
    inner_join(top_MeanRatio_genes |>
                   filter(MeanRatio > 1) |>
                  select(cell_type = cellType.target, gene)) |>
    group_by(cell_type) |>
    summarise(overlap_top10 = paste(gene, collapse = ", ")))

write_csv(marker_summary, file = here(data_dir, sprintf("subtype_marker_summary_%s.csv", celltype)))

#### Enrichment dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_enrichment_genes$test[match(rownames(sce), top_enrichment_genes$gene)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_enrichment.pdf", celltype)))
sce |>
    scDotPlot(features = top_enrichment_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "Marker" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

#### spatial registration ####
# 
# layer_modeling_results <- list()
# 
# ## ERC
# layer_modeling_results$spatialERC <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))
# colnames(layer_modeling_results$spatialERC$enrichment) <- gsub("_Sp", "~Sp", colnames(layer_modeling_results$spatialERC$enrichment))
# 
# ## PsychENCODE enrichment
# layer_modeling_results$snDLPFC_PEC <- list()
# layer_modeling_results$snDLPFC_PEC$enrichment <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/registration_stats_LIBD.rds")
# 
# ## sestan sn enrichment
# layer_modeling_results$sestan_EC <- readRDS(here("processed-data", "04_snRNA-seq", "24_external_data_check", "sestan_EC_modeling.rds"))
# 
# ## spatial HPC modeling 
# layer_modeling_results$spatialHPC <- list()
# layer_modeling_results$spatialHPC$enrichment <- read.csv("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)
# 
# names(layer_modeling_results)
# # [1] "HumanPilot"   "spatialDLPFC" "spatialERC"   "snDLPFC_PEC"  "sestan_EC"
# 
# cor_layer <- map(layer_modeling_results, function(layer_mod){
#     cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
#                                 modeling_results = layer_mod,
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

#### Oligo Grubman Cor ####
if(celltype == "Oligo"){
    
    grubman_subtypes_fn <- list.files(here("external-data", "Grubman2019"), pattern = "Oligo_", full.names = TRUE)
    names(grubman_subtypes_fn) <- gsub("Grubman_Oligo_|.csv", "", basename(grubman_subtypes_fn))
    
    grubman_subtypes_data <- map2_dfr(grubman_subtypes_fn, names(grubman_subtypes_fn), ~read.csv(.x) |> 
                                          mutate(g_cell_type = .y)) |>
        rename(gene = geneName)
    
    grubman_subtypes_data |> count(g_cell_type)
    
    enrichment_genes <- sig_genes_extract(
        n = nrow(sce_pseudo),
        modeling_results = modeling_results,
        model_type = "enrichment",
        reverse = FALSE,
        sce_layer = sce_pseudo,
        gene_name = "gene_name"
    )

    erc_v_grubman <- enrichment_genes |>
        inner_join(grubman_subtypes_data, relationship = "many-to-many")
    
    erc_v_grubman_cor <- erc_v_grubman |>
        group_by(test, g_cell_type) |>
        summarise(n = n(),
                  cor = cor(logFC, LFC))
    
    write_csv(erc_v_grubman_cor, file = here(data_dir, "erc_v_grubman_oligo_cor.csv"))
    
    erc_v_grubman_cor |>
        group_by(test) |> 
        arrange(-cor) |>
        slice(1)
    
    # test    g_cell_type     n   cor
    # <chr>   <chr>       <int> <dbl>
    # 1 Oligo.1 o5             95 0.706
    # 2 Oligo.2 o5             95 0.544
    # 3 Oligo.3 o4            651 0.752
    # 4 Oligo.4 o5             95 0.358
    # 5 Oligo.5 o4            651 0.482
    
    erc_v_grubman_cor |>
        group_by( g_cell_type) |> 
        arrange(-cor) |>
        slice(1) 
    
    # test    g_cell_type     n   cor
    # <chr>   <chr>       <int> <dbl>
    # 1 Oligo.3 o1            411 0.483
    # 2 Oligo.2 o2            394 0.397
    # 3 Oligo.1 o3             96 0.303
    # 4 Oligo.3 o4            651 0.752
    # 5 Oligo.1 o5             95 0.706
    # 6 Oligo.1 o6            212 0.511
    
    (erc_v_grubman_cor_wide <- erc_v_grubman_cor |>
            select(-n) |>
            pivot_wider(names_from = "test", values_from = "cor") |>
            column_to_rownames("g_cell_type") |>
            as.matrix())
    
    # Oligo.1    Oligo.2    Oligo.3     Oligo.4    Oligo.5
    # o1 -0.4213429 -0.3937171  0.4830757 -0.08971312  0.2299773
    # o2  0.3170668  0.3972655 -0.4272965  0.03675172 -0.2132549
    # o3  0.3034454  0.2883027 -0.4272468  0.28193107 -0.2031632
    # o4 -0.6556795 -0.6273667  0.7524666 -0.33494686  0.4816174
    # o5  0.7063523  0.5439999 -0.6512348  0.35844107 -0.6110251
    # o6  0.5113262  0.4239644 -0.5416082  0.12368402 -0.3155310
    
    
    ## grubman annotations 
    grubmab_anno <- data.frame(annotation = c("AD", "AD","AD", "Undetermined", "Control", "Control"))
    rownames(grubmab_anno) <- paste0("o", 1:6)
    
    grubman_col_ha <- HeatmapAnnotation(df = grubmab_anno,
                                         col = list("annotation" = c(AD = "purple",Undetermined = "black", Control = "Green")))
    
    pdf(here(plot_dir, "Oligo_grubman_cor_logFC.pdf"), height = 4, width = 8)
    Heatmap(t(erc_v_grubman_cor_wide), name = "logFC cor",
            bottom_annotation = grubman_col_ha)
    dev.off()
        
    
    ## grubman marker dot plot
    
    grubman_o_makers <- grubman_subtypes_data |>
        group_by(g_cell_type) |>
        arrange(-LFC) |>
        slice(1:5) |>
        ungroup() |>
        group_by(gene) |>
        slice_max(LFC) |>
        filter(gene %in% rownames(sce)) |>
        arrange(g_cell_type)
    
    # grubman_o_makers |>
    #     ungroup() |>
    #     count(gene) |>
    #     filter(n==2)
    #     
    
    rowData(sce)$Grubman_Oligo <- NULL
    rowData(sce)$Grubman_Oligo <- grubman_o_makers$g_cell_type[match(rownames(sce), grubman_o_makers$gene)] 
    table(rowData(sce)$Grubman_Oligo)
    
    pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_Grubman_Oligo.pdf", celltype)))
    sce |>
        scDotPlot(features = grubman_o_makers$gene,
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "Grubman_Oligo",
                  scale = TRUE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno),
                  clusterRows = FALSE,
                  groupLegends = FALSE)
    dev.off()
    

#### Oligo Sadick Cor ####
    
    sadick_stab3 <- readxl::read_xlsx(here("external-data", "Sadick2022","Sadick2022_STab3.xlsx"), sheet = "LEN_so_oligo_DEGs") |>
        mutate(Sadick_cluster = paste0("Sadick_S", LEN_so_oligo_cluster),
               data = "Sadick")
    
    sadick_stab5 <- readxl::read_xlsx(here("external-data", "Sadick2022","Sadick2022_STab5.xlsx"), sheet = "Oligo_merged_DEGs") |>
        mutate(Sadick_cluster = paste0("Sadick_Int", Oligo_merged_cluster),
               data = "Integrated")
    
    sadick_Oligo_markers <- sadick_stab3 |> bind_rows(sadick_stab5)
    
    sadick_Oligo_markers |> dplyr::count(Sadick_cluster)
    # Sadick_cluster     n
    # <chr>          <int>
    # 1 Sadick_Ol0         4
    # 2 Sadick_Ol1        88
    # 3 Sadick_Ol2        37
    # 4 Sadick_Ol3        95
    # 5 Sadick_Ol4        53
    # 6 Sadick_Ol5       332
    # 7 Sadick_Ol6       250
    
    erc_v_sadick <- enrichment_genes |>
        inner_join(sadick_Oligo_markers, relationship = "many-to-many")
    
    erc_v_sadick_cor <- erc_v_sadick |>
        group_by(test, data, Sadick_cluster) |>
        summarise(n = n(),
                  cor = cor(logFC, avg_log2FC))
    
    write_csv(erc_v_sadick_cor, file = here(data_dir, "erc_v_sadick_oligo_cor.csv"))
    
    erc_v_sadick_cor |>
        group_by(test, data) |> 
        slice_max(cor)
    
    erc_v_sadick_cor |>
        group_by(Sadick_cluster, data) |> 
        slice_max(cor)
     
    (erc_v_sadick_cor_wide <- erc_v_sadick_cor |>
            select(-n) |>
            group_by(data) |>
            group_map(~pivot_wider(.x, names_from = "test", values_from = "cor") |>
                          column_to_rownames("Sadick_cluster") |>
                          as.matrix())
    )
    
    ## heatmap
    pdf(here(plot_dir, "Oligo_Sadick_cor_logFC.pdf"), height = 4, width = 8)
    map(erc_v_sadick_cor_wide, ~print(Heatmap(t(.x), name = "logFC cor")))
    dev.off()
        
    
    ## sadick marker dot plot
    map2(list(sadick_stab3, sadick_stab5), c("S", "Int"), function(data, name){
        
        sadick_o_makers <- data |>
            group_by(Sadick_cluster) |>
            arrange(-avg_log2FC) |>
            slice(1:5) |>
            filter(gene %in% rownames(sce)) 
        
        # sadick_o_makers |>
        #     ungroup() |>
        #     count(gene) |>
        #     filter(n==2)
        # #     
        # 
        # sadick_o_makers |> filter(gene %in% c("RASGRF1", "RBFOX1"))
        
        rowData(sce)$sadick_Oligo <- NULL
        rowData(sce)$sadick_Oligo <- sadick_o_makers$Sadick_cluster[match(rownames(sce), sadick_o_makers$gene)] 
        table(rowData(sce)$sadick_Oligo)
        
        pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_sadick%s_Oligo.pdf", celltype, name)))
        print(sce |>
            scDotPlot(features = unique(sadick_o_makers$gene),
                      group = "cell_type_anno",
                      groupAnno = "cell_type_anno",
                      featureAnno = "sadick_Oligo",
                      scale = TRUE,
                      annoColors = list("cell_type_anno" = cell_type_colors$anno),
                      clusterRows = FALSE,
                      groupLegends = FALSE))
        
        print(sce |>
            scDotPlot(features = unique(sadick_o_makers$gene),
                      group = "cell_type_anno",
                      groupAnno = "cell_type_anno",
                      featureAnno = "sadick_Oligo",
                      scale = TRUE,
                      annoColors = list("cell_type_anno" = cell_type_colors$anno),
                      clusterRows = TRUE,
                      groupLegends = FALSE))
        
        dev.off()
    })
    
    #### Oligo Siletti2023 Cor ####
    # Table S5: DEGs between Oligo1 (OPALIN) and Oligo2 (RBFOX).
    Siletti2023 <- read_csv(here("external-data", "Siletti2023", "science.add7046_table_s5.csv")) |>
        rename(gene = `...1`, Siletti_stat = stat)
    
    Siletti2023_cor <- Siletti2023 |>
        inner_join(enrichment_genes |> 
                       filter(top <= 100),
                   relationship = "many-to-many") |>
        group_by(test) |> 
        summarise(n = n(),
                  cor_t_stat = cor(stat, Siletti_stat),
                  cor_logFC = cor(log2FoldChange, logFC))
    
    erc_v_Siletti_cor_wide <- Siletti2023_cor |>
        select(-n) |>
        column_to_rownames("test") |>
                      as.matrix()
    
    pdf(here(plot_dir, "Oligo_Siletti2023_cor_logFC.pdf"), height = 4, width = 5)
    print(Heatmap(erc_v_Siletti_cor_wide, 
                  name = "cor", 
                  cluster_columns = FALSE,
                  column_title = "Oligo1 (OPALIN) vs. Oligo2 (RBFOX)"))
    dev.off()
    
    
}


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

# disease associated 
#### Cell type checks ####

enrichment_genes <- sig_genes_extract(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)


if(celltype == "Oligo"){
    #### Oligo ####
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
        pdf_fn = here(plot_dir, "sn_violin_lit_markers.pdf"),
        color_pal = cell_type_colors$anno
    )
    
    ## Bernie genes of interest
    oligo_bernie_genes <- plot_gene_express(
        sce,
        genes = c("PSEN1", "NCSTN"),
        category = "cell_type_anno",
        color_pal = cell_type_colors$anno
    )
    
    ggsave(oligo_bernie_genes, filename = here(plot_dir, "sn_violin_Oligo_check.png"))
    
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
                  annoColors = list("cell_type_anno" = cell_type_colors$anno),
                  clusterRows = FALSE,
                  groupLegends = FALSE)
    dev.off()
    
    #### Oligo Marques markers ####
    
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
        color_pal = cell_type_colors$anno
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
                  annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                    "Marques_markers" = Marques_palette),
                  clusterRows = FALSE,
                  groupLegends = FALSE)
    
    sce |>
        scDotPlot(features = unique(Marques_markers),
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "Marques_markers",
                  scale = TRUE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                    "Marques_markers" = Marques_palette),
                  clusterRows = TRUE,
                  groupLegends = FALSE)
    dev.off()
    
    
} else if(celltype == "Astro"){
    #### Astro ####
    lit_markers <- list(disease_associated = c("SERPINA3", "C4B", "TNFRSF1A", "IL1B", "IL33", "HMOX1", "TNF", "ERK", "ERK2"), #https://doi.org/10.1038/s41593-025-01873-x
                        AD_risk = c("APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1"),
                        Mathys_subtypes = c("GRM3","DPP10", "DCLK1" ,"LUZP2")
    )
    
    
    lit_markers <- map(lit_markers, ~.x[.x %in% rownames(sce)])
    
    plot_marker_express_List(
        sce,
        gene_list = lit_markers,
        cellType_col = "cell_type_anno",
        pdf_fn = here(plot_dir, "sn_violin_lit_markers.pdf"),
        color_pal = cell_type_colors$anno
    )
    

    ## lit gene dot plot
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
                  scale = FALSE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno),
                  clusterRows = FALSE,
                  groupLegends = FALSE)
    dev.off()
    
    astro_subtype_markers <- list(
        protoplasmic = c("SLC1A2", "SLC1A3", "GLUL", "GPC4", "NDRG2"),
        fibrous      = c("GFAP", "VIM", "ID3", "ID4", "HOPX"),
        perivascular = c("CP", "AQP4", "KCNJ10", "AGRN", "CD44")
    )
    
    
    astro_subtype_markers <- map(astro_subtype_markers, ~.x[.x %in% rownames(sce)])
    
    plot_marker_express_List(
        sce,
        gene_list = astro_subtype_markers,
        cellType_col = "cell_type_anno",
        pdf_fn = here(plot_dir, "sn_violin_Astro_suntype_markers.pdf"),
        color_pal = cell_type_colors$anno
    )
    

    ## lit gene dot plot
    astro_subtype_markers <- AnnotationDbi::unlist2(astro_subtype_markers)
    
    rowData(sce)$AstroMarker <- NULL
    rowData(sce)$AstroMarker <- names(astro_subtype_markers)[match(rownames(sce), astro_subtype_markers)] 
    table(rowData(sce)$AstroMarker)
    
    pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_subtype.pdf", celltype)))
    sce |>
        scDotPlot(features = astro_subtype_markers,
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "AstroMarker",
                  scale = TRUE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno),
                  clusterRows = FALSE,
                  groupLegends = TRUE)
    dev.off()

     #### Astro cor Green 2024 ####       
    
    Green_Astro_stats <- readxl::read_xlsx(here("external-data", "Green2024","Green2024_SuppTable2.xlsx"), sheet = "DEGs") |>
        filter(cell.type == "astrocytes",
        !is.na(avg_log2FC)) |>
        mutate(Green_Astro = state)

    Green_Astro_stats |> count(cell.type)
        Green_Astro_stats |> count(Green_Astro)

    erc_v_Green <- enrichment_genes |>
        inner_join(Green_Astro_stats, relationship = "many-to-many")

    erc_v_Green_cor <- erc_v_Green |>
        group_by(test, Green_Astro) |>
        summarise(n = n(),
                cor = cor(logFC, avg_log2FC))

    write_csv(erc_v_Green_cor, file = here(data_dir, "erc_v_Green2024_Astro_cor.csv"))

    erc_v_Green_cor |>
        group_by(test) |> 
        slice_max(cor)

    erc_v_Green_cor |>
        group_by(Green_Astro) |> 
        slice_max(cor)


    (erc_v_Green_cor_wide <- erc_v_Green_cor |>
            select(-n) |>
            pivot_wider(names_from = "test", values_from = "cor") |>
            column_to_rownames("Green_Astro") |>
            as.matrix())

    pdf(here(plot_dir, "Astro_Green_cor_logFC.pdf"), height = 4, width = 8)
    Heatmap(t(erc_v_Green_cor_wide), name = "logFC cor")
    dev.off()
    
    #### Astro cor Mathys2024 ####       
    Mathys_Astro_stats <- readxl::read_xlsx(here("external-data", "Mathys2024", "41586_2024_7606_MOESM3_ESM","Supplementary_Table_2_Gene_expression_marker_lists.xlsx"), sheet = "Marker_genes_Ast") |>
        mutate(Mathys_Astro = gsub(" ", "_", cluster)) |>
        select(-Column1)

    Mathys_Astro_stats |> count(Mathys_Astro)

    erc_v_Mathys <- enrichment_genes |>
        inner_join(Mathys_Astro_stats, relationship = "many-to-many")

    erc_v_Mathys_cor <- erc_v_Mathys |>
        group_by(test, Mathys_Astro) |>
        summarise(n = n(),
                cor = cor(logFC, avg_log2FC))

    write_csv(erc_v_Mathys_cor, file = here(data_dir, "erc_v_Mathys2024_Astro_cor.csv"))

    erc_v_Mathys_cor |>
        group_by(test) |> 
        slice_max(cor)

    erc_v_Mathys_cor |>
        group_by(Mathys_Astro) |> 
        slice_max(cor)


    (erc_v_Mathys_cor_wide <- erc_v_Mathys_cor |>
            select(-n) |>
            pivot_wider(names_from = "test", values_from = "cor") |>
            column_to_rownames("Mathys_Astro") |>
            as.matrix())

    pdf(here(plot_dir, "Astro_Mathys_cor_logFC.pdf"), height = 4, width = 8)
    Heatmap(t(erc_v_Mathys_cor_wide), name = "logFC cor")
    dev.off()
    
}else if(celltype == "Micro"){
    #### Micro ####
    ## from https://www.biocompare.com/Editorial-Articles/590587-A-Guide-to-Oligodendrocyte-Markers/
    lit_markers <- list(disease_associated = c("TREM2"), #https://doi.org/10.1038/s41593-025-01873-x
                        inflamation = c("NLPR3")) 
    
    lit_markers <- map(lit_markers, ~.x[.x %in% rownames(sce)])
    
    plot_marker_express_List(
        sce,
        gene_list = lit_markers,
        cellType_col = "cell_type_anno",
        pdf_fn = here(plot_dir, sprintf("sn_violin_%s_markers.pdf", celltype)),
        color_pal = cell_type_colors$anno
    )
    
    # ## oligo marker dot plot
    # lit_markers <- AnnotationDbi::unlist2(lit_markers)
    # 
    # rowData(sce)$Marker <- NULL
    # rowData(sce)$Marker <- names(lit_markers)[match(rownames(sce), lit_markers)] 
    # table(rowData(sce)$Marker)
    # 
    # pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_lit.pdf", celltype)))
    # sce |>
    #     scDotPlot(features = lit_markers,
    #               group = "cell_type_anno",
    #               groupAnno = "cell_type_anno",
    #               featureAnno = "Marker",
    #               scale = TRUE,
    #               annoColors = list("cell_type_anno" = cell_type_colors$anno),
    #               clusterRows = FALSE,
    #               groupLegends = FALSE)
    # dev.off()
    # 
} else if(celltype == "OPC") {
  #### OPC  Marques markers ####  
    Marques_markers_OPC <- list(OPC = c("PDGFRA", "CSPG4", "PTPRZ1", "PCDH15"), 
                             COP = c("VCAN",  "SOX6", "GPR17", "NEU4", "BMP4", "NKX2-2"),
                             Siletti = c("HES5", "IRX5", "GPC5"))
    
    Marques_markers_OPC <- map(Marques_markers_OPC, ~.x[.x %in% rownames(sce)])
    
    Marques_markers_OPC <- AnnotationDbi::unlist2(Marques_markers_OPC)
    
    rowData(sce)$Marques_markers_OPC <- NULL
    rowData(sce)$Marques_markers_OPC <- names(Marques_markers_OPC)[match(rownames(sce), Marques_markers_OPC)] 
    table(rowData(sce)$Marques_markers_OPC)
    
    Marques_palette_OPC <- c(
        OPC   = "navy",
        COP  = "pink")
    
    
    pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_Marques_markers_OPC.pdf", celltype)))
    sce |>
        scDotPlot(features = unique(Marques_markers_OPC),
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "Marques_markers_OPC",
                  scale = TRUE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                    "Marques_markers_OPC" = Marques_palette_OPC),
                  clusterRows = FALSE,
                  groupLegends = FALSE)
    
    sce |>
        scDotPlot(features = unique(Marques_markers_OPC),
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "Marques_markers_OPC",
                  scale = TRUE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                    "Marques_markers_OPC" = Marques_palette_OPC),
                  clusterRows = TRUE,
                  groupLegends = FALSE)

    dev.off()
    
}


# slurmjobs::job_single('35_sn_subcluster_marker_modeling', create_shell = TRUE, memory = '25G', command = "Rscript 35_sn_subcluster_marker_modeling.R --celltype Oligo")

slurmjobs::job_loop(loops = list(celltype = c("Astro", "Micro", "Endo", "OPC", "Vasc", "Excit", "Inhib")), create_shell = TRUE, name = "35_sn_subcluster_marker_modeling", create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

