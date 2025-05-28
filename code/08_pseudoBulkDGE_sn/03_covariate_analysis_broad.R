## Louise Huuki-Myers, April 2025
## Examine covairates in sn dataset for DGE from broad cell type

library("SingleCellExperiment")
library("scran")
library("tidyverse")
library("ggrepel")
library("ggpubr")
# library("variancePartition")
library("here")
library("sessioninfo")

#### Set up dirs ####
# data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "03_covariate_analysis_broad")
#if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Set up dirs ####
plot_dir <- here("plots", "08_pseudoBulkDGE_sn", "03_covariate_analysis_broad")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
# message(Sys.time(), " - Load HDF5 sce")
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS"))

## load colors
cell_type_colors <- metadata(sce_pb)$cell_type_colors
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

colnames(colData(sce_pb))

table(sce_pb$cell_type_broad)
# Astro Macro Micro   OPC Oligo  Vasc Excit Inhib 
# 31    19    31    31    31    30    31    31 

# Extract the relevant columns from the data
pd <- as.data.frame(colData(sce_pb)) |>
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE,
           cell_type_broad, ncells, exp_round, seq_round, Rin,
           pseudo_sum_umi, pseudo_expr_chrM, pseudo_expr_chrM_ratio) 


#### Generate n sample barplots for each cell_type ####

cell_type_count <- pd |> group_by(APOE_carrier, cell_type_broad) |> count()

pb_cell_type_bar_APOE_carrier <- pd |>
    count(APOE_carrier, cell_type_broad) |> 
    ggplot(aes(x = cell_type_broad, y=n, fill = APOE_carrier)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()

ggsave(pb_cell_type_bar_APOE_carrier, filename = here(plot_dir, 'pb_cell_type_broad_bar_APOE_carrier.png'))

pb_cell_type_bar_APOE <- pd |>
    count(APOE, cell_type_broad) |> 
    ggplot(aes(x = cell_type_broad, y=n, fill = APOE)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()

ggsave(pb_cell_type_bar_APOE, filename = here(plot_dir, 'pb_cell_type_broad_bar_APOE.png'))

#### t-test variables ####

test_variables <- c("ncells", "pseudo_sum_umi", "pseudo_expr_chrM_ratio")
names(test_variables) <- test_variables

summary(pd[,test_variables])

var_t_test <- map(test_variables, function(test_var){
    
    y_position <- pd |>
        group_by(cell_type_broad) |>
        summarise(y.position = max(!!sym(test_var)) + .05*max(!!sym(test_var)))
    
    var_t_test <- pd |>
        do(compare_means(!!sym(test_var) ~ APOE_carrier, data = ., method = "t.test", p.adjust.method = "fdr", group.by = "cell_type_broad")) |>
        ungroup() |>
        mutate(p.signif.fdr = case_when(p.adj < 0.005 ~ "***",
                                        p.adj < 0.01 ~"**",
                                        p.adj < 0.05 ~"*",
                                        TRUE~""),
               fdr_anno = sprintf("FDR=%.3f%s", p.adj, p.signif.fdr)) |>
        left_join(y_position)
    
    
    
    boxplot_test_var <- pd |>
        ggplot() +
        geom_boxplot(aes(y = !!sym(test_var), x = APOE_carrier, fill = APOE_carrier), outlier.shape = NA) +
        geom_point(aes(y = !!sym(test_var), x = APOE_carrier, fill = APOE_carrier)) +
        geom_text(data = var_t_test, aes(label = fdr_anno, x = 1, y = y.position))+
        facet_wrap(~cell_type_broad, scales = "free_y") +
        scale_fill_manual(values = APOE_carrier_colors) +
        theme_bw() 
    
    ggsave(boxplot_test_var, filename = here(plot_dir, sprintf("boxplot_broad_%s_APOE_carrier.png", test_var)), width = 10)
    
    return(var_t_test)
    
})

map(var_t_test, ~.x|> filter(p.adj < 0.05))


#### variance partition ####

# load sample level data
dim(sce_pb)
table(sce_pb$APOE)
# E2/E2 E2/E3 E3/E4 E4/E4 
# 46    62    73    54 
table(sce_pb$APOE_carrier)
# E2+ E4+ 
# 108 127

any(is.na(pd$Rin))

#### Variance Partition data ####
varPart_summary <- map_dfr(list.files(here("processed-data", "08_pseudoBulkDGE_sn", "02_VariancePartition_sn"), full.names = TRUE),
                           read.csv) 

any(is.na(varPart_summary$cell_type))

head(varPart_summary)

varPart_summary |> count(cell_type)
varPart_summary |> filter(metric == "Median", name != "Residuals") |> group_by(cell_type, form) |> slice_max(value)
varPart_summary |> filter(metric == "Mean", name != "Residuals")  |> group_by(cell_type, form) |> slice_max(value)

varPart_heatmap <- function(my_form = "carrier", my_metric = "Mean"){
    
    varPart_summary_matrix <- varPart_summary |> 
        filter(metric == my_metric, form == my_form,  name != "Residuals") |>
        select(cell_type, name, value) |>
        pivot_wider(values_from = "value", names_from = "cell_type") |>
        column_to_rownames("name") |>
        as.matrix()
    
    ## order by rowSums    
    vp_metric_sum <- rowSums(varPart_summary_matrix)
    varPart_summary_matrix <- varPart_summary_matrix[order(vp_metric_sum, decreasing = TRUE), ]
    
    col_fun = colorRamp2(c(0, max(varPart_summary_matrix)), c("white", "red"))
    
    pdf(here(plot_dir, sprintf("sn_varPart_summary-%s-%s.pdf", my_form, my_metric)), height = 8, width = 11)
    print(Heatmap(varPart_summary_matrix,
                  name = my_metric,
                  col = col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = TRUE,
                  column_title = sprintf("Visium VarPart - %s", my_form)
    ))
    dev.off()
    
    rownames(varPart_summary_matrix) <- paste0(my_form, "-", rownames(varPart_summary_matrix))
    
    return(varPart_summary_matrix)
}

## plot heatmaps
(form_list <- unique(varPart_summary$form))
names(form_list) <- form_list

var_part_mean <- map(form_list, ~varPart_heatmap(my_form = .x, my_metric = "Mean"))

## make combined heatmap
var_part_mean_all <- do.call("rbind", var_part_mean)

covar_index <- order(rowSums(var_part_mean_all))

var_part_mean_all <- var_part_mean_all[covar_index,]

model_row_ha_simple <- rowAnnotation(
    mod = jaffelab::ss(rownames(var_part_mean_all), "-"),
    covar = jaffelab::ss(rownames(var_part_mean_all), "-",2)
)

col_fun = colorRamp2(c(0, max(var_part_mean_all)), c("white", "red"))
pdf(here(plot_dir, "sn_varPart_summary-ALL-Mean.pdf"), height = 12, width = 11)
print(Heatmap(var_part_mean_all,
              name = "Mean",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              right_annotation = model_row_ha_simple,
              column_title = sprintf("sn VarPart - %s", "ALL")
))
dev.off()


## plot medians
walk(form_list, ~varPart_heatmap(my_form = .x, my_metric = "Median"))

#### colinearity ####

log_fn <- list.files(here("code", "08_pseudoBulkDGE_sn", "logs"), pattern = "02_VariancePartition_sn", full.names = TRUE)
names(log_fn) <- gsub(".txt", "",gsub("02_VariancePartition_Visium_", "", basename(log_fn)))
logs <- map(log_fn, readLines)

colinear_pairs <- map(logs, ~.x[(grep("High colinearity between variables:", .x)+1):(grep("null device ", .x)-1)])
map_int(colinear_pairs, length)

all(map_lgl(colinear_pairs[-1], ~setequal(.x, colinear_pairs[[1]])))

colinear_pairs[[1]]

# slurmjobs::job_single('03_covariate_analysis_broad', create_shell = TRUE, memory = '50G', command = "Rscript 03_covariate_analysis_broad.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()