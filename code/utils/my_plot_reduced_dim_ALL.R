
my_plot_reduced_dim_ALL <- function(prefix, suffix, nested_dir = TRUE){
    
    if(nested_dir){
        plot_dir <- here(plot_dir,paste0("reduced_dims_", suffix))
        if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    } 
    
    source(here("code", "utils", "my_plot_reduced_dim.R"))
    ## swap rownames to gene names for plotting markers
    rownames(sce) <- rowData(sce)$Symbol 
    
    ## categorical
    walk(rdim_name = c("UMAP", "TSNE"){
        walk(c("sample_id", "seq_round", "exp_round","APOE","quick_cluster"), ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", sufix = suffix))
        
        ## continuous
        walk(c("sum", "detected", "subsets_Mito_percent"), ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "con", sufix = suffix))
        
        ## Expression 
        walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(sce, prefix = preix, dimred = rdim_name, var_type = "express", my_var = .x, sufix = suffix))
    })

    
}