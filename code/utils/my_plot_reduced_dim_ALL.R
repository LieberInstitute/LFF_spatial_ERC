source(here("code", "utils", "my_plot_reduced_dim.R"))
load(here::here("processed-data", "project_colors.Rdata"))

my_plot_reduced_dim_ALL <- function(prefix, suffix, nested_dir = TRUE){
    
    if(nested_dir){
        plot_dir <- here(plot_dir,paste0("reduced_dims_", suffix))
        if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    } 
    
    ## swap rownames to gene names for plotting markers
    rownames(sce) <- rowData(sce)$Symbol 
    
    ## categorical
    walk(c("UMAP", "TSNE"), function(rdim_name){
        
        ## Categorical vars 
        cat_vars <- c("seq_round", "exp_round", "quick_cluster")
        if(any(!cat_vars %in% colnames(colData(sce)))) warning("Missing cat vars: ", paste0(cat_vars[!cat_vars %in% colnames(colData(sce))], collapse = ", "))
        
        cat_vars <- cat_vars[cat_vars %in% colnames(colData(sce))]
        
        walk(cat_vars, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix))
        walk(cat_vars, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix, facet = TRUE))
        
        
        ## Categorical vars with colors
        cat_vars <- c("Ancestry","sex","APOE","quick_cluster", "sample_id")
        cat_colors <- c(ancestry_colors, sex_colors, APOE_genotype_colors, sample_colors)
        
        if(any(!cat_vars %in% colnames(colData(sce)))) warning("Missing cat vars: ", paste0(cat_vars[!cat_vars %in% colnames(colData(sce))], collapse = ", "))
        cat_vars <- cat_vars[cat_vars %in% colnames(colData(sce))]
        cat_colors <- cat_colors[cat_vars %in% colnames(colData(sce))]
        
        walk2(cat_vars, cat_colors, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix, color_pal = .y))
        walk2(cat_vars, cat_colors, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix, facet = TRUE, color_pal = .y))
        
        ## Continuous variables
        con_vars <- c("sum", "detected", "subsets_Mito_percent")
        if(any(!con_vars %in% colnames(colData(sce)))) warning("Missing con vars: ", paste0(con_vars[!con_vars %in% colnames(colData(sce))], collapse = ", "))
        con_vars <- con_vars[con_vars %in% colnames(colData(sce))]
        
        walk(con_vars,  ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "con", suffix = suffix))
        
        ## Expression 
        genes <- c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119")
        if(any(!genes %in% rownames(rowData(sce)))) warning("Missing genes: ", paste0(genes[!genes %in% rownames(rowData(sce))], collapse = ", "))
        genes <- genes[genes %in% rownames(rowData(sce))]
        
        walk(genes, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, var_type = "express", my_var = .x, suffix = suffix))
    })
}