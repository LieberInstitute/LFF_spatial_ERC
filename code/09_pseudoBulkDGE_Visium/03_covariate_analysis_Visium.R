## Louise Huuki-Myers, April 2025
## Examine covaraites in Visium dataset for DGE
## adapted from https://github.com/LieberInstitute/dlpfc_asd/blob/a250a1d7e20bd754c5f1186aa96ce0752d55e556/code/08_pseudoBulkDGE_s/02_covariate_analysis.R

library("spatialLIBD")
library("SingleCellExperiment")
library("tidyverse")
library("ggpubr")
library("ggrepel")
library("GGally")
library("here")
library("sessioninfo")
# library("variancePartition")
library("ComplexHeatmap")
library("circlize")

#### Set up dirs ####
plot_dir <- here("plots", "09_pseudoBulkDGE_Visium", "03_covariate_analysis_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
spe_pb <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))

## load colors
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

# Extract the relevant columns from the data
pd <- as.data.frame(colData(spe_pb))

table(spe_pb$SpD)

spd_tab <- pd |> count(SpD_syn, SpD) |> select(-n)

#### t-test variables ####

test_variables <- c("ncells", "pseudo_sum_umi", "pseudo_expr_chrM_ratio")
names(test_variables) <- test_variables

summary(pd[,test_variables])

var_t_test <- map(test_variables, function(test_var){
    
    y_position <- pd |>
        group_by(SpD) |>
        summarise(y.position = max(!!sym(test_var)) + .05*max(!!sym(test_var)))
    
    var_t_test <- pd |>
        do(compare_means(!!sym(test_var) ~ APOE_carrier, data = ., method = "t.test", p.adjust.method = "fdr", group.by = "SpD")) |>
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
        facet_wrap(~SpD, scales = "free_y") +
        scale_fill_manual(values = APOE_carrier_colors) +
        theme_bw() 
    
    ggsave(boxplot_test_var, filename = here(plot_dir, sprintf("Visium_boxplot_%s_APOE_carrier.png", test_var)), width = 10)
    
    return(var_t_test)
    
})

map(var_t_test, ~.x|> filter(p.adj < 0.05))

#### Explore quality metrics ####

scater_umi_v_chrM <- pd |> 
    ggplot(aes(x = pseudo_sum_umi, y = pseudo_expr_chrM, color = SpD)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = SpD_colors)
    
ggsave(scater_umi_v_chrM, filename = here(plot_dir, "spe_pb_scater_umi_v_chrM.png"))

scater_umi_v_chrM_ratio <- pd |> 
    ggplot(aes(x = pseudo_sum_umi, y = pseudo_expr_chrM_ratio, color = SpD)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = SpD_colors)
    
ggsave(scater_umi_v_chrM_ratio, filename = here(plot_dir, "spe_pb_scater_umi_v_chrM_ratio.png"))

scater_ncells_v_umi <- pd |> 
    ggplot(aes(x = ncells, y = pseudo_sum_umi, color = SpD)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = SpD_colors)
    
ggsave(scater_ncells_v_umi, filename = here(plot_dir, "spe_pb_scater_ncells_v_umi.png"))

scater_ncells_v_chrM_ratio <- pd |> 
    ggplot(aes(x = ncells, y = pseudo_expr_chrM_ratio, color = SpD)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = SpD_colors)
    
ggsave(scater_ncells_v_chrM_ratio, filename = here(plot_dir, "spe_pb_scater_ncells_v_chrM_ratio.png"))

pd |> 
    group_by(SpD) |>
    summarise(cor_cell_umi = cor(ncells, pseudo_sum_umi),
              cor_cell_chrM = cor(ncells, pseudo_expr_chrM_ratio),
              cor_umi_chrM = cor(pseudo_sum_umi, pseudo_expr_chrM_ratio))

ggpair_plot <- ggpairs(pd, 
                       columns = c("ncells", "pseudo_sum_umi","pseudo_expr_chrM","pseudo_expr_chrM_ratio"), 
                       aes(colour = SpD)) +
    scale_color_manual(values = SpD_colors) +
    scale_fill_manual(values = SpD_colors) +
    theme_bw()

ggsave(ggpair_plot, filename = here(plot_dir, "spe_pb_ggpair_quality_metrics.png"), height = 10, width = 10)



#### Variance Partition data ####
varPart_summary <- map_dfr(list.files(here("processed-data", "09_pseudoBulkDGE_Visium", "02_VariancePartition_Visium"), full.names = TRUE),
                              read.csv) |>
    rename(SpD_syn = SpD) |> 
    left_join(spd_tab) 

any(is.na(varPart_summary$SpD))

head(varPart_summary)

varPart_summary |> count(SpD)

varPart_summary |> filter(metric == "Median", name != "Residuals") |> group_by(SpD, form) |> slice_max(value)
varPart_summary |> filter(metric == "Mean", name != "Residuals") |> head()


varPart_heatmap <- function(my_form = "carrier", my_metric = "Mean"){
    
    varPart_summary_matrix <- varPart_summary |> 
        filter(metric == my_metric, form == my_form,  name != "Residuals") |>
        select(SpD, name, value) |>
        pivot_wider(values_from = "value", names_from = "SpD") |>
        column_to_rownames("name") |>
        as.matrix()
    
        col_fun = colorRamp2(c(0, max(varPart_summary_matrix)), c("white", "red"))
    
    pdf(here(plot_dir, sprintf("Visium_varPart_summary-%s-%s.pdf", my_form, my_metric)), height = 8, width = 11)
    print(Heatmap(varPart_summary_matrix,
            name = my_metric,
            col = col_fun,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            column_title = sprintf("Visium VarPart - %s", my_form)
    ))
    dev.off()
    
    return(varPart_summary_matrix)
    }

## plot heatmaps
(form_list <- unique(varPart_summary$form))

var_part_mean <- map(form_list, ~varPart_heatmap(my_form = .x, my_metric = "Mean"))

walk(form_list, ~varPart_heatmap(my_form = .x, my_metric = "Median"))

#### colinearity ####

log_fn <- list.files(here("code", "09_pseudoBulkDGE_Visium", "logs"), pattern = "02_VariancePartition_Visium", full.names = TRUE)
names(log_fn) <- gsub(".txt", "",gsub("02_VariancePartition_Visium_", "", basename(log_fn)))
logs <- map(log_fn, readLines)

colinear_pairs <- map(logs, ~.x[(grep("High colinearity between variables:", .x)+1):(grep("null device ", .x)-1)])
map_int(colinear_pairs, length)

all(map_lgl(colinear_pairs[-1], ~setequal(.x, colinear_pairs[[1]])))

colinear_pairs[[1]]
# [1] "  APOE and APOE_num"                    "  APOE and sample_id"                  
# [3] "  APOE_num and sample_id"               "  sample_id and Sex"                   
# [5] "  sample_id and Age"                    "  sample_id and Anc_Afr"               
# [7] "  sample_id and Rin"                    "  sample_id and Visium_slide"          
# [9] "  sample_id and round"                  "  sample_id and ncells"                
# [11] "  sample_id and pseudo_sum_umi"         "  sample_id and pseudo_expr_chrM_ratio"
# [13] "  Visium_slide and round"    

# slurmjobs::job_single('03_covariate_analysis_Visium', create_shell = TRUE, memory = '10G', command = "Rscript 03_covariate_analysis_Visium.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()