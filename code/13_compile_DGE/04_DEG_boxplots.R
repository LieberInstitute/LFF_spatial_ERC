## Louise Huuki-Myers, July 2025
## Plot boxplots for DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("DeconvoBuddies")
library("SingleCellExperiment")
library("patchwork")

source(here("code","utils","plot_DEG_express.R"))
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype <- "sn_fine"
# opt$datatype <- "sn_broad"
# opt$datatype <- "Visium"

print(opt)

#### set up dirs ####
# data_dir <- here("processed-data", "13_compile_DGE", "04_DEG_boxplots")
# if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "04_DEG_boxplots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load pseudobulk data ####
if(opt$datatype == "Visium"){
    pb_fn <- here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS")
    cluster_var <- "SpD"
} else if(opt$datatype == "sn_broad"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS")
    cluster_var <- "cell_type_broad"
}else if(opt$datatype == "sn_fine"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS")
    cluster_var <- "cell_type_anno"
} else {
    stop("non-valid datatype")
}

message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, basename(pb_fn)))
sce_pb <- readRDS(pb_fn)
dim(sce_pb)
table(sce_pb$registration_variable)

## Drop sample Br1289
sce_pb <- sce_pb[,sce_pb$BrNum != 'Br1289']

rownames(sce_pb) <- rowData(sce_pb)$gene_name

load(here("processed-data", "project_colors.Rdata"))

#### Carrier DE data ####

DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
DE_data <- readRDS(DE_data_fn)

DEGs_signif <- DE_data |> 
    group_by(cluster) |>
    arrange(vlmf_adj.P.Val) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 ~ "up",
                                vlmf_logFC < 0 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class),
           rank = row_number())

DEGs_signif |> filter(vlmf_adj.P.Val < 0.05) |> dplyr::count(cluster)

## check model
head(model.matrix(~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio, colData(sce_pb)))

cluster_levels <- levels(sce_pb[[cluster_var]])
cluster_levels2 <- as.character(unique(DE_data$cluster))

cluster_levels <- intersect(cluster_levels, cluster_levels2)

all(cluster_levels %in% sce_pb[[cluster_var]])

# plot_DEG_express(sce = sce_pb,
#                     stats = DE_data,
#                     clus = "Astro",
#                     n_genes = 10,
#                     pval_col = "vlmf_adj.P.Val",
#                     fc_col = "vlmf_logFC",
#                     gene_col = "gene_name",
#                     cluster_col = cluster_var,
#                     category_col = "APOE_carrier",
#                     mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
#                     color_pal = APOE_carrier_colors,
#                     plot_points = FALSE,
#                     ncol = 2,
#                     cleanY_P = 4)

table(sce_pb[[cluster_var]])

pdf(here(plot_dir, sprintf("DEG_boxplots_carrier_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express(sce = sce_pb,
                                          stats = DE_data,
                                          clus = .x,
                                          n_genes = 10,
                                          pval_col = "vlmf_adj.P.Val",
                                          fc_col = "vlmf_logFC",
                                          gene_col = "gene_name",
                                          cluster_col = cluster_var,
                                          category_col = "APOE_carrier",
                                          mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                          color_pal = APOE_carrier_colors,
                                          plot_points = TRUE,
                                          ncol = 2,
                                          cleanY_P = 4)
         )
dev.off()

pdf(here(plot_dir, sprintf("DEG_boxplots_carrier_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express(sce = sce_pb,
                                      stats = DE_data,
                                      clus = .x,
                                      n_genes = 10,
                                      pval_col = "vlmf_adj.P.Val",
                                      fc_col = "vlmf_logFC",
                                      gene_col = "gene_name",
                                      cluster_col = cluster_var,
                                      category_col = "APOE_carrier",
                                      mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                      color_pal = APOE_carrier_colors,
                                      plot_points = TRUE,
                                      ncol = 2,
                                      cleanY_P = 4)
)
dev.off()

if(opt$datatype == "sn_broad"){

    # DE_data |> filter(cluster == "Oligo", gene_name %in% c("PLP1", "VCAM1"))
    # DE_data |> filter(cluster == "Oligo") |> arrange(vlmf_adj.P.Val)
    # DE_data |> filter(cluster == "Oligo", vlmf_adj.P.Val < 0.05) |> arrange(-vlmf_logFC)
        
    #Oligo
    map(c("OPALIN", "SRARP", "CD69","C1QL3", "NPTXR"), function(g){
        dge_plot <- plot_DEG_express(
                sce = sce_pb,
                stats = DE_data,
                gene = g,
                cluster_col = "cell_type_broad",
                clus = "Oligo",
                gene_col = "gene_name",
                color_pal = APOE_carrier_colors,
                plot_points = TRUE
        ) +
            labs(title = "Oligo")
        ggsave(dge_plot, 
               filename = here(plot_dir, sprintf("DEG_boxplots_carrier_%s_%s_%s.png", opt$datatype, "Oligo", g)),
               height = 3, width = 3)
        
    })
    
    #Astro
    map(c("C1QL3", "NPTXR"), function(g){
        dge_plot <- plot_DEG_express(
                sce = sce_pb,
                stats = DE_data,
                gene = g,
                cluster_col = "cell_type_broad",
                clus = "Oligo",
                gene_col = "gene_name",
                color_pal = APOE_carrier_colors,
                plot_points = TRUE
        ) +
            labs(title = "Oligo")
        ggsave(dge_plot, 
               filename = here(plot_dir, sprintf("DEG_boxplots_carrier_%s_%s_%s.png", opt$datatype, "Oligo", g)),
               height = 4, width = 3)
        
    })

}else if(opt$datatype == "sn_fine"){
    
    #Oligo.3
    map(c("MAPT", "FOS"), function(g){
        dge_plot <- plot_DEG_express(
                sce = sce_pb,
                stats = DE_data,
                gene = g,
                cluster_col = "cell_type_anno",
                clus = "Oligo.3",
                gene_col = "gene_name",
                color_pal = APOE_carrier_colors,
                plot_points = TRUE
        ) +
            labs(title = "Oligo.3")
        ggsave(dge_plot, 
               filename = here(plot_dir, sprintf("DEG_boxplots_carrier_%s_%s_%s.png", opt$datatype, "Oligo.3", g)),
               height = 4, width = 3)
        
    })

} else if(opt$datatype == "Visium"){
    
    #Oligo
    boxplots_klk6 <- map(c("Vasc~Sp09D08", "WM.uf~Sp09D07"), function(spd){
        dge_plot <- plot_DEG_express(
            sce = sce_pb,
            stats = DE_data,
            gene = "KLK6",
            cluster_col = "SpD",
            clus = spd,
            gene_col = "gene_name",
            color_pal = APOE_carrier_colors,
            plot_points = TRUE
        ) +
            labs(title = spd)
        return(dge_plot)
    })
    ggsave(boxplots_klk6[[1]] | boxplots_klk6[[2]],
           filename = here(plot_dir, sprintf("DEG_boxplots_carrier_%s_%s.png", opt$datatype, "KLK6")),
           height = 3, width = 6)
}


#### plot taupath data ####

## load tau pathology data
tau_tb <- read.csv(here("processed-data", "00_project_prep","05_pathology","sample_taupathy.csv")) |>
    column_to_rownames("BrNum")

sce_pb$taupathy <- ifelse(tau_tb[sce_pb$BrNum,]$taupathy, "t+", "t-")
table(sce_pb$taupathy)

sce_pb$carrier_tau <- paste(sce_pb$APOE_carrier, sce_pb$taupathy)

table(sce_pb$registration_variable, sce_pb$carrier_tau)

carrier_tau_colors <- c(`E2+ t-` = "#398A84",
                        `E2+ t+` = "#60BEB8",
                        `E4+ t-` = "#D46B43",
                        `E4+ t+` = "#DD8A69")

pdf(here(plot_dir, sprintf("DEG_boxplots_carrier_taupathy_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express(sce = sce_pb,
                                         stats = DE_data,
                                         clus = .x,
                                         n_genes = 10,
                                         pval_col = "vlmf_adj.P.Val",
                                         fc_col = "vlmf_logFC",
                                         gene_col = "gene_name",
                                         cluster_col = cluster_var,
                                         category_col = "carrier_tau",
                                         mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                         color_pal = carrier_tau_colors,
                                         plot_points = TRUE,
                                         ncol = 2,
                                         cleanY_P = 4)
)
dev.off()

#### Carrier by Ancestry DE data ####

DE_ancestry_data_fn <- here("processed-data", "13_compile_DGE",  "05_compile_DGE_ancestry", opt$datatype, sprintf("DGE_results_ancestry_%s.Rds", opt$datatype))
DE_ancestry_data <- readRDS(DE_ancestry_data_fn)

## update cluster levels
cluster_levels <- levels(sce_pb[[cluster_var]])
cluster_levels2 <- as.character(unique(DE_ancestry_data$cluster))
cluster_levels <- intersect(cluster_levels, cluster_levels2)


DE_ancestry_data |> filter(vlmf_adj.P.Val < 0.05) |> dplyr::count(contrast, cluster)

sce_pb$carrier_Anc <- factor(paste(sce_pb$APOE_carrier, sce_pb$Ancestry), levels = c("E2+ AA","E4+ AA", "E2+ EA", "E4+ EA"))
table(sce_pb$carrier_Anc)

ancestry_mod <- model.matrix(~0+carrier_Anc + Sex + Age + pseudo_expr_chrM_ratio, colData(sce_pb))
colnames(ancestry_mod)


pdf(here(plot_dir, sprintf("DEG_boxplots_ancestry_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express_contrast_top(sce = sce_pb,
                                               stats = DE_ancestry_data,
                                               clus = .x,
                                               n_genes = 10,
                                               pval_col = "vlmf_adj.P.Val",
                                               fc_col = "vlmf_logFC",
                                               gene_col = "gene_name",
                                               cluster_col = cluster_var,
                                               category_col = "carrier_Anc",
                                               contrast_1 = "carrier_AA",
                                               contrast_2 = "carrier_EA",
                                               mod = ancestry_mod,
                                               plot_points = TRUE,
                                               ncol = 2,
                                               cleanY_P = 4)
)
dev.off()

anc_single_gene <- function(gene_list, 
                            cluster_col = registration_variable){
    
    map2(gene_list, names(gene_list), function(genes, cl){
        
        map(genes, function(g){
            dge_plot <- plot_DEG_express_contrast(
                sce = sce_pb,
                stats = DE_ancestry_data,
                gene = g,
                cluster_col = cluster_col,
                clus = cl,
                gene_col = "gene_name",
                category_col = "carrier_Anc",
                # color_pal = APOE_carrier_colors,
                plot_points = TRUE
            ) 
            ggsave(dge_plot, 
                   filename = here(plot_dir, sprintf("DEG_boxplots_ancestry_%s_%s_%s.png", opt$datatype, cl, g)),
                   height = 3, width = 3)
        })
    })
}



if(opt$datatype == "Visium"){
    anc_single_gene(gene_list = list(`WM_Sp09D06` = c("PEX14")))
    
} else if(opt$datatype == "sn_broad"){
    
    anc_single_gene(gene_list = list(Astro = c("TCTN1", "JUP"),
                                     Vasc = c("RIPK2")), 
                    cluster_col = "cell_type_broad")
}


#### Carrier by Sex DE data ####

DE_Sex_data_fn <- here("processed-data", "13_compile_DGE",  "09_compile_DGE_Sex", opt$datatype, sprintf("DGE_results_Sex_%s.Rds", opt$datatype))
DE_Sex_data <- readRDS(DE_Sex_data_fn)

cluster_levels <- levels(sce_pb[[cluster_var]])
cluster_levels2 <- as.character(unique(DE_Sex_data$cluster))
cluster_levels <- intersect(cluster_levels, cluster_levels2)

DE_Sex_data |> filter(vlmf_adj.P.Val < 0.05) |> dplyr::count(contrast, cluster)

sce_pb$carrier_Sex <- factor(paste(sce_pb$APOE_carrier, sce_pb$Sex), levels = c("E2+ F","E4+ F", "E2+ M", "E4+ M"))
table(sce_pb$carrier_Sex)

Sex_mod <- model.matrix(~0+carrier_Sex + Ancestry + Age + pseudo_expr_chrM_ratio, colData(sce_pb))
colnames(Sex_mod)

DE_Sex_data |> 
    filter(vlmf_adj.P.Val < 0.05) |>
    dplyr::count(contrast, cluster)

DE_Sex_data |> 
    group_by(cluster, contrast) |>
    arrange(vlmf_adj.P.Val) |>
    dplyr::slice(1) |>
    select(cluster, contrast, gene_name, vlmf_adj.P.Val, vlmf_logFC)

pdf(here(plot_dir, sprintf("DEG_boxplots_Sex_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express_contrast_top(sce = sce_pb,
                                                   stats = DE_Sex_data,
                                                   clus = .x,
                                                   n_genes = 10,
                                                   pval_col = "vlmf_adj.P.Val",
                                                   fc_col = "vlmf_logFC",
                                                   gene_col = "gene_name",
                                                   cluster_col = cluster_var,
                                                   category_col = "carrier_Sex",
                                                   contrast_1 = "carrier_F",
                                                   contrast_2 = "carrier_M",
                                                   mod = Sex_mod,
                                                   plot_points = TRUE,
                                                   ncol = 2,
                                                   cleanY_P = 4)
)
dev.off()

sex_single_gene <- function(gene_list, 
                            cluster_col = registration_variable){
    
    map2(gene_list, names(gene_list), function(genes, cl){
        
        map(genes, function(g){
            dge_plot <- plot_DEG_express_contrast(
                sce = sce_pb,
                stats = DE_Sex_data,
                gene = g,
                cluster_col = cluster_col,
                clus = cl,
                gene_col = "gene_name",
                category_col = "carrier_Sex",
                contrast_1 = "carrier_F",
                contrast_2 = "carrier_M",
                # color_pal = APOE_carrier_colors,
                plot_points = TRUE
            ) 
            ggsave(dge_plot, 
                   filename = here(plot_dir, sprintf("DEG_boxplots_Sex_%s_%s_%s.png", opt$datatype, cl, g)),
                   height = 3, width = 3)
        })
    })
}



if(opt$datatype == "Visium"){
    anc_single_gene(gene_list = list(`WM_Sp09D06` = c("PEX14")))
    
} else if(opt$datatype == "sn_broad"){
    
    sex_single_gene(gene_list = list(Oligo = c("ADAMTS10", "RIOK2"),
                                     Vasc = c("CCN3"),
                                     Macro = c("CCL4L2", "EZR")), 
                    cluster_col = "cell_type_broad")
}


#### Interaction DE data ####

DE_interaction_data_fn <- here("processed-data", "13_compile_DGE",  "05_compile_DGE_interaction", opt$datatype, sprintf("DGE_results_interaction_%s.Rds", opt$datatype))
DE_interaction_data <- readRDS(DE_interaction_data_fn)

DE_interaction_data |> filter(vlmf_adj.P.Val < 0.2) |> dplyr::count(cluster)
DE_interaction_data |> filter(vlmf_adj.P.Val < 0.05) 

table(sce_pb$carrier_Anc)

interaction_mod <- model.matrix(~APOE_carrier_syn*Ancestry + Sex + Age + pseudo_expr_chrM_ratio, colData(sce_pb))
interaction_mod <-interaction_mod[,c("(Intercept)", "APOE_carrier_synE4:AncestryEA", "APOE_carrier_synE4","AncestryEA", "SexM","Age","pseudo_expr_chrM_ratio")]
head(interaction_mod)

## resubset clusters
cluster_levels2 <- as.character(unique(DE_interaction_data$cluster))
cluster_levels <- intersect(cluster_levels, cluster_levels2)


pdf(here(plot_dir, sprintf("DEG_boxplots_interaction_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express(sce = sce_pb,
                                      stats = DE_interaction_data,
                                      clus = .x,
                                      n_genes = 10,
                                      pval_col = "vlmf_adj.P.Val",
                                      fc_col = "vlmf_logFC",
                                      gene_col = "gene_name",
                                      cluster_col = cluster_var,
                                      category_col = "carrier_Anc",
                                      mod = interaction_mod,
                                      # color_pal = carrier_tau_colors,
                                      plot_points = TRUE,
                                      ncol = 2,
                                      cleanY_P = 2)
)
dev.off()

# slurmjobs::job_single('04_DEG_boxplots', create_shell = TRUE, memory = '5G', command = "Rscript 04_DEG_boxplots.R --datatype sn_broad")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

