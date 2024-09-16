## September 2024, Louise Huuki-Myers
## Compare ERC spatial Domains to DLPFC layers

library("spatialLIBD")
library("purrr")
library("here")
library("sessioninfo")

#### Set up dirs ####
plot_dir <- here("plots", "05_spe_correct_cluster", "10_spatial_registration_DLPFC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "10_spatial_registration_DLPFC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### get reference layer enrichment statistics ####
layer_modeling_results <- map(c(HumanPilot = "modeling_results", spatialDLPFC = "spatialDLPFC_Visium_modeling_results"), fetch_data)

#### load data ####
modeling_fn <- list.files(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_PCA_Harmony"),
           pattern = "modeling_results",
           full.names = TRUE)

names(modeling_fn) <- gsub("modeling_results-|.rds", "", basename(modeling_fn))

cor_humanPilot <- map2(modeling_fn, names(modeling_fn), function(fn, name){
    
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
        
        pdf(here(plot_dir, sprintf("layer_stat_cor-%s_vs_ERC-%s.pdf", dataset, name)))
        layer_stat_cor_plot(cor_layer, max = max(cor_layer))
        dev.off()
        
        return(cor_layer)
        
    })
    
    
    anno <- map(cor_layer, ~annotate_registered_clusters(
        cor_stats_layer = .x,
        confidence_threshold = 0.25,
        cutoff_merge_ratio = 0.25
    ))
    
    return(list(cor_layer = cor_layer, layer_anno = anno))
})


## spatialDLPFC spatial domain annotations 
"processed-data/rdata/spe/08_spatial_registration/bayesSpace_layer_annotations.csv"

