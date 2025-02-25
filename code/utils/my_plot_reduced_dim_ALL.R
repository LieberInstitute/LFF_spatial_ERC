source(here("code", "utils", "my_plot_reduced_dim.R"))
load(here::here("processed-data", "project_colors.Rdata"))

my_plot_reduced_dim_ALL <- function(prefix, suffix, nested_dir = TRUE, plot_dir_rd){
    
    if(nested_dir){
        plot_dir_rd <- here(plot_dir,paste0("reduced_dims_", suffix))
        if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    } 
    
    message("Plot dir: ", plot_dir_rd)
    
    ## swap rownames to gene names for plotting markers
    rownames(sce) <- rowData(sce)$Symbol 
    
    ## Categorical vars 
    cat_vars <- c("seq_round", "exp_round", "quick_cluster")
    if(any(!cat_vars %in% colnames(colData(sce)))) message("Missing cat vars: ", paste0(cat_vars[!cat_vars %in% colnames(colData(sce))], collapse = ", "))
    
    cat_vars <- cat_vars[cat_vars %in% colnames(colData(sce))]
    
    ## Categorical vars with colors
    cat_vars_color <- c("Ancestry","Sex","APOE", "sample_id")
    cat_colors <- list(Ancestry = ancestry_colors, 
                    Sex = sex_colors, 
                    APOE = APOE_genotype_colors,
                    sample_id = sample_colors)
    
    if(any(!cat_vars_color %in% colnames(colData(sce)))) message("Missing color cat vars: ", paste0(cat_vars_color[!cat_vars_color %in% colnames(colData(sce))], collapse = ", "))
    cat_vars_color <- cat_vars_color[cat_vars_color %in% colnames(colData(sce))]
    cat_colors <- cat_colors[cat_vars_color %in% colnames(colData(sce))]
    
    ## Continuous variables
    con_vars <- c("sum", "detected", "subsets_Mito_percent")
    if(any(!con_vars %in% colnames(colData(sce)))) message("Missing con vars: ", paste0(con_vars[!con_vars %in% colnames(colData(sce))], collapse = ", "))
    con_vars <- con_vars[con_vars %in% colnames(colData(sce))]
    
    ## Expression 
    genes <- c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119")
    if(any(!genes %in% rownames(rowData(sce)))) message("Missing genes: ", paste0(genes[!genes %in% rownames(rowData(sce))], collapse = ", "))
    genes <- genes[genes %in% rownames(rowData(sce))]
    
    ## Plot
    walk(c("UMAP", "TSNE"), function(rdim_name){
        
        ## cat vars
        message("** Plot cat variables - ", rdim_name)
        walk(cat_vars, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix, plot_dir = plot_dir_rd))
        walk(cat_vars, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix, facet = TRUE, plot_dir = plot_dir_rd))
        
        ## cat vars colors
        message("** Plot cat variables w/ color - ", rdim_name)
        walk2(cat_vars_color, cat_colors, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix, color_pal = .y, plot_dir = plot_dir_rd))
        walk2(cat_vars_color, cat_colors, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "cat", suffix = suffix, facet = TRUE, color_pal = .y, plot_dir = plot_dir_rd))
        
        ## continuous
        message("** Plot continuous variables - ", rdim_name)
        walk(con_vars,  ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, my_var = .x, var_type = "con", suffix = suffix, plot_dir = plot_dir_rd))
        
        ## gene expression
        message("** Plot gene expression - ", rdim_name)
        walk(genes, ~my_plot_reduced_dim(sce, prefix = prefix, dimred = rdim_name, var_type = "express", my_var = .x, suffix = suffix, plot_dir = plot_dir_rd))
    })
}
