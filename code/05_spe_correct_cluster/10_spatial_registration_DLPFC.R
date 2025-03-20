## September 2024, Louise Huuki-Myers
## Compare ERC spatial Domains to DLPFC layers

library("spatialLIBD")
library("purrr")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "10_spatial_registration_DLPFC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### get reference layer enrichment statistics ####
layer_modeling_results <- map(c(HumanPilot = "modeling_results", spatialDLPFC = "spatialDLPFC_Visium_modeling_results"), fetch_data)

#### load data ####
modeling_fn <- c(
    list.files(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm"),
           pattern = "modeling_results",
           full.names = TRUE),
    list.files(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_Markers"),
           pattern = "modeling_results",
           full.names = TRUE)
    )

names(modeling_fn) <- gsub("modeling_results-|.rds", "", basename(modeling_fn))

modeling_fn

## Add spatialDLPFC spatial domain annotations to modeling
dlpfc_anno <- read.csv("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_spatial_registration/bayesSpace_layer_annotations.csv") |>
    dplyr::filter(bayesSpace == "k09") |>
    mutate(layer_combo2 = gsub(" ", "~", layer_combo2))

dlpfc_colnames <- colnames(layer_modeling_results$spatialDLPFC$enrichment)
pwalk(dlpfc_anno, function(...) dlpfc_colnames <<- gsub(..5, ..3, dlpfc_colnames))
colnames(layer_modeling_results$spatialDLPFC$enrichment) <- dlpfc_colnames
head(layer_modeling_results$spatialDLPFC$enrichment)


cor_anno <- map2(modeling_fn, names(modeling_fn), function(fn, name){
    
    message(Sys.time()," - ", name)
    
    erc_modeling_results <- readRDS(fn)
    
    ## extract t-statics and rename
    registration_t_stats <- erc_modeling_results$enrichment[, grep("^t_stat", colnames(erc_modeling_results$enrichment))]
    colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))
    
    
    cor_layer <- map2(layer_modeling_results, names(layer_modeling_results), function(layer_mod, dataset){
        
        cor_layer <- layer_stat_cor(stats = registration_t_stats,
                                    modeling_results = layer_mod,
                                    model_type = "enrichment",
                                    top_n = 100)
        
        cor_layer <- cor_layer[,order(colnames(cor_layer))]
        
        ## plot
        # pdf(here(plot_dir, sprintf("layer_stat_cor-%s_vs_ERC-%s.pdf", dataset, name)))
        # layer_stat_cor_plot(cor_layer, max = max(cor_layer))
        # dev.off()
        
        return(cor_layer)
        
    })
    
    anno <- map(cor_layer, ~annotate_registered_clusters(
        cor_stats_layer = .x,
        confidence_threshold = 0.25,
        cutoff_merge_ratio = 0.25
    ))
    
    return(list(cor_layer = cor_layer, layer_anno = anno))
})

## save data
save(cor_anno, file = here(data_dir, "spatial_registration_erc_v_DLPFC_cor_anno.Rdata"))

cor_anno$BayesSpace_SVGm_k09$cor$spatialDLPFC

cor_anno$BayesSpace_SVGm_k09$layer_anno$HumanPilot
# cluster layer_confidence   layer_label layer_label_simple
# 1 Sp09D06             good            WM                 WM
# 2 Sp09D07             good            WM                 WM
# 3 Sp09D05             good        Layer1                 L1
# 4 Sp09D08             good        Layer1                 L1
# 5 Sp09D04             good        Layer6                 L6
# 6 Sp09D03             good Layer5/Layer4               L5/4
# 7 Sp09D09             good        Layer4                 L4
# 8 Sp09D01             good Layer2/Layer3               L2/3
# 9 Sp09D02             good        Layer3                 L3

cor_anno$BayesSpace_SVGm_k09$layer_anno$spatialDLPFC
# cluster layer_confidence                                 layer_label
# 1 Sp09D06             good                                  WM~Sp09D06
# 2 Sp09D07             good                       WM~Sp09D06/WM~Sp09D09
# 3 Sp09D05             good                                  L1~Sp09D02
# 4 Sp09D08             good                                  L1~Sp09D01
# 5 Sp09D04             good                                  L6~Sp09D07
# 6 Sp09D03             good L5~Sp09D04/L4~Sp09D08/L3~Sp09D05/L6~Sp09D07
# 7 Sp09D02             good                       L3~Sp09D05/L2~Sp09D03
# 8 Sp09D01             good                       L3~Sp09D05/L2~Sp09D03
# 9 Sp09D09             good L4~Sp09D08/L3~Sp09D05/L5~Sp09D04/L2~Sp09D03

#### complex heatmap ####
# 
# erc_colors <- c(Sp09D01 = "#F6222E", #Vivid_Red
#                    Sp09D02 = "#FE00FA", #Vivid_Purple
#                    Sp09D03 = "#16FF32", #Vivid_Yellowish_Green
#                    Sp09D04 = "#3283FE", #Strong_Purplish_Blue
#                    Sp09D05 = "#FEAF16", #Vivid_Orange_Yellow
#                    Sp09D06 = "#E4E1E3", #Purplish_White
#                    Sp09D07 = "#5A5156", #Dark_Purplish_Gray
#                    Sp09D08 = "#B00068", #Vivid_Purplish_Red
#                    Sp09D09 = "#1CFFCE") #Brilliant_Green
# 
# erc_colors <- c(ERC_D1 = "#F6222E", #Vivid_Red
#                 ERC_D2 = "#FE00FA", #Vivid_Purple
#                 ERC_D3 = "#16FF32", #Vivid_Yellowish_Green
#                 ERC_D4 = "#3283FE", #Strong_Purplish_Blue 
#                 ERC_D5 = "#FEAF16", #Vivid_Orange_Yellow
#                 ERC_D6 = "#E4E1E3", #Purplish_White
#                 ERC_D7 = "#5A5156", #Dark_Purplish_Gray
#                 ERC_D8 = "#B00068", #Vivid_Purplish_Red
#                 ERC_D9 = "#1CFFCE") #Brilliant_Green


walk2(cor_anno, names(cor_anno), function(cor, name){
    pdf(here(plot_dir, sprintf("spatial_registration_heatmap_erc_%s_v_%s.pdf", knice, dataset)), height = 5, width = 5)
    print(layer_stat_cor_plot_complex(
        cor_stat_layer = cor$cor_layer,
        annotation = cor$layer_anno
        
    ))
    dev.off()
})
# })
    
    
walk2(cor_anno, names(cor_anno), function(ca, name){

    #get numeric k
    k <- readr::parse_number(name)
    knice <- gsub("BayesSpace_SVGm_", "", name)
        message(Sys.time(), " - ", knice)
    
    ## create color pallet
    erc_colors <- Polychrome::palette36.colors(k)
    if(k ==2) erc_colors <- erc_colors[seq(2)] ## fix return 3 colors bug
    
    names(erc_colors) <- sort(rownames(ca$cor_layer[[1]]))
    
    pwalk(list(cor_matrix = ca$cor_layer, 
               dataset=names(ca$cor_layer), 
               anno= ca$layer_anno), 
          function(cor_matrix, dataset, anno){
        
        ## spatial domain color bar
        erc_colors <- erc_colors[rownames(cor_matrix)]
        
        if(dataset == "HumanPilot"){
            ref_color <- spatialLIBD::libd_layer_colors
            } else ref_color <- NULL
        
        # Create heatmap
        pdf(here(plot_dir, sprintf("spatial_registration_heatmap_erc_%s_v_%s.pdf", knice, dataset)), height = 5, width = 5)
        print(layer_stat_cor_plot_complex(
            cor_stats_layer = cor_matrix,
            query_colors = erc_colors,
            reference_colors = ref_color,
            annotation = anno,
            cluster_columns = FALSE
            
        ))
        dev.off()
    })
    
})

# slurmjobs::job_single('10_spatial_registration_DLPFC', create_shell = TRUE, memory = '25G', command = "Rcript 10_spatial_registration_DLPFC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()