## Louise Huuki-Myers, June 2025
## Compile and plot dreamlet sn data

library("SpatialExperiment")
library("dreamlet")
library("tidyverse")
library("here")
library("sessioninfo")
library("DFplyr")
library("ggrepel")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "05_compile_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "05_compile_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

#### Load data ####
## processed dreamlet input
res.proc <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))

## dreamlet results
res.dl.fn <- list.files(here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn"), full.names = TRUE)
names(res.dl.fn) <- gsub(".RDS", "", gsub("dreamlet_sn-", "", basename(res.dl.fn)))

res.dl.list <- purrr::map(res.dl.fn, readRDS)
# purrr::map(res.dl.list, coefNames)

coef_list <- list(apoe_i = c("APOEE2/E3:Anc_Afr", "APOEE3/E4:Anc_Afr", "APOEE4/E4:Anc_Afr"),
                  apoe_n0 = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  apoe = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  carrier_i = "APOE_carrierE4+:Anc_Afr",
                  carrier_n0 = "APOE_carrierE4+",
                  carrier_kr = "APOE_carrierE4+",
                  carrier_n1 = "APOE_carrierE4+",
                  carrier_n2 = "APOE_carrierE4+",
                  carrier_n3 = "APOE_carrierE4+",
                  carrier_n4 = "APOE_carrierE4+",
                  carrier_sf = "APOE_carrierE4+",
                  carrier = "APOE_carrierE4+",
                  e4e4_i = "APOE_E4E4TRUE:Anc_Afr",
                  e4e4_n0 = "APOE_E4E4TRUE",
                  e4e4 = "APOE_E4E4TRUE")

coef_list_kr <- coef_list
names(coef_list_kr) <- paste0(names(coef_list_kr), "_kr")

coef_list <- c(coef_list, coef_list_kr)

## check all coef are present, order correctly
all(names(res.dl.list) %in% names(coef_list))
coef_list <- coef_list[names(res.dl.list)]

#### TopTable ####
## get top table
tt <- map2(res.dl.list, coef_list, ~topTable(.x, .y, number = Inf))

map_int(tt, nrow)

## correct pvalue by cluster, add model
tt <- map2(tt, names(tt), ~.x |> 
               group_by(assay) |>
               mutate(adj.P.Val.cell_type = p.adjust(P.Value),
                      model = .y))

map(tt, ~.x |> summarise(global_FDR10 = sum(adj.P.Val < 0.1),
                         global_FDR20 = sum(adj.P.Val < 0.2),
                         cell_type_FDR10 = sum(adj.P.Val.cell_type < 0.1),
                         cell_type_FDR20 = sum(adj.P.Val.cell_type < 0.2)
                         )
    )

main_mods <- c("carrier", "apoe", "e4e4", "carrier_i", "apoe_i", "e4e4_i")

fdr_model_count <- map2_dfr(tt, names(tt), ~.x |> 
                                ungroup() |>
                                summarise(global_FDR10 = sum(adj.P.Val < 0.1),
                                          global_FDR20 = sum(adj.P.Val < 0.2),
                                          cell_type_FDR10 = sum(adj.P.Val.cell_type < 0.1),
                                          cell_type_FDR20 = sum(adj.P.Val.cell_type < 0.2)) |>
                                mutate(model = .y) |>
                                as.data.frame())  |>
    mutate(main_mod = model %in% main_mods) |>
    arrange(-main_mod)

write.csv(fdr_model_count, file = here(data_dir, "dreamlet_sn_FDR_model_count.csv"), row.names = FALSE)

## filter and export key results
tt_signif <- map_dfr(tt[main_mods], ~.x |> as.data.frame() |> filter(adj.P.Val.cell_type < 0.2))

tt_signif |> select(assay, ID, model) |> group_by(model, assay) |> summarise(n = n(), DEGs = paste(ID, collapse = ", "))

write.csv(tt_signif, file = here(data_dir, "dreamlet_sn_topTable_FDR20.csv"), row.names = FALSE)

## save main mods
tt_sn <- tt[main_mods]

map_int(tt_sn, ~.x |> filter(adj.P.Val.cell_type < 0.2))

save(tt_sn, file = here(data_dir, "dreamlet_sn_TopTables.Rdata"))
# load(here("processed-data", "10_dreamlet_sn", "05_compile_dreamlet_sn", "dreamlet_sn_TopTables.Rdata"))

#### bar plots ###
# pval_bar <- tt_sn[["apoe"]] |>
#     mutate(adj.P.Val.bin = cut(adj.P.Val.cell_type, breaks = 0:10/10)) |>
#     ggplot(aes(x = assay, fill = adj.P.Val.bin)) +
#     geom_bar()
# 
# ggsave(pval_bar, filename = here(plot_dir, "pval_bar_test.png"))

#### check ddr models ####
fdr_model_count |>
    mutate(ddr = ifelse(grepl("kr", model), "kr", "s"),
           model = gsub("_kr", "", model)) |>
    arrange(model) |>
    group_by(model) |>
    filter("kr" %in% ddr)

# global_FDR10 cell_type_FDR10 model      ddr  
# <int>           <int> <chr>      <chr>
# 1            0               1 apoe       kr   
# 2            0               4 apoe       s    
# 3            0               0 carrier    kr   
# 4            0               1 carrier    s    
# 5            0               0 carrier_sf kr   
# 6            0               1 carrier_sf s    
# 7            0               0 e4e4       kr   
# 8            0               1 e4e4       s  

tt_kr_comapre <- tt[["carrier"]] |>
    as.data.frame() |>
    bind_rows(as.data.frame(tt[["carrier_kr"]])) |>
    select(assay, ID, t, model) |>
    pivot_wider(values_from = t, names_from = model) |>
    mutate(diff = carrier - carrier_kr)

summary(tt_kr_comapre$carrier - tt_kr_comapre$carrier_kr)

tt_kr_comapre_scatter <- tt_kr_comapre |>
    ggplot(aes(carrier, carrier_kr, color = assay)) +
    geom_point(size = 1, alpha = .5) +
    geom_abline() +
    theme_bw() +
    labs(title = "carrier model t-stats")

ggsave(tt_kr_comapre_scatter, filename = here(plot_dir, "tt_kr_comapre_scatter.png"))


#### pval distribution ####
walk(names(tt), function(mod){
    
    pval_histo <- tt[[mod]] |>
        ggplot(aes(x = P.Value)) +
        geom_histogram(binwidth = 0.01) +
        facet_wrap(~assay) +
        theme_bw() +
        labs(title = mod)
    
    ggsave(pval_histo, filename = here(plot_dir, sprintf("dreamlet_sn_pval_histo-%s.png", mod)), width = 10)
    
})

walk(main_mods, function(mod){
    
    pval_histo <- tt[[mod]] |>
        ggplot(aes(x = adj.P.Val)) +
        geom_histogram(binwidth = 0.01) +
        facet_wrap(~assay) +
        theme_bw() +
        labs(title = mod)
    
    ggsave(pval_histo, filename = here(plot_dir, sprintf("dreamlet_sn_pval_histo_adj-%s.png", mod)), width = 10)
    
})

walk(main_mods, function(mod){
    
    pval_histo <- tt[[mod]] |>
        ggplot(aes(x = adj.P.Val, y = adj.P.Val.cell_type)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_text_repel(aes(label = ifelse(adj.P.Val < 0.3 | adj.P.Val.cell_type < 0.3, ID, "")), size = 1.5) +
        facet_wrap(~assay) +
        theme_bw() +
        geom_hline(yintercept = 0.2, linetype ="dashed", color = "blue") +
        geom_vline(xintercept = 0.2, linetype ="dashed", color = "red") +
        labs(title = mod)
    
    ggsave(pval_histo, filename = here(plot_dir, sprintf("dreamlet_sn_pval_adj_scatter-%s.png", mod)), width = 10)
    
})


#### check sex DGE ####
tt_sex <- topTable(res.dl.list[["carrier_sf"]], "SexM", number = Inf)

tt_sex |>
    group_by(assay) |>
    mutate(adj.P.Val.cell_type = p.adjust(P.Value)) |>
    summarise(global_FDR10 = sum(adj.P.Val < 0.1),
              cell_type_FDR10 = sum(adj.P.Val.cell_type < 0.1))

# DataFrame with 8 rows and 3 columns
# groups[xx, -ncol(groups)] global_FDR10 cell_type_FDR10
# <character>    <integer>       <integer>
# 1                     Astro           25              10
# 2                     Excit           18              10
# 3                     Inhib           13               5
# 4                     Macro            1               1
# 5                     Micro            7               3
# 6                     Oligo           16               6
# 7                       OPC           12               4
# 8                      Vasc            7               5

## check VarPart genes
tt[["carrier"]] |>
    as.data.frame() |>
    filter(ID %in% c("CLOCK", "ABCC4"))

tt[["apoe"]] |>
    as.data.frame() |>
    filter(ID %in% c("ADAM10", "ADAM17"))

tt[["e4e4"]] |>
    as.data.frame() |>
    filter(ID %in% c("SRGAP2C", "NAALADL2", "APP"))

#### Volcano Plots ####
pdf(here(plot_dir, "sn_dreamlet_VolcanoPlot.pdf"), height = 11, width = 8)
# plotVolcano(res.dl, coef = "APOE_carrierE4+")

map(c("carrier","e4e4","carrier_i","e4e4_i"), ~plotVolcano(res.dl.list[[.x]], coef_list[[.x]], cutoff = 0.2) + labs(title = .x))

dev.off()

# p_limit <- tt[["carrier"]] |>
#     as.data.frame() |>
#     filter(adj.P.Val.cell_type < 0.2) |>
#     group_by(assay) |>
#     summarise(p_limit = max(P.Value))
# 
# 
# tt[["carrier"]] |>
#         mutate(DE_class = case_when(adj.P.Val.cell_type < 0.2 ~ "E4+",
#                                     adj.P.Val.cell_type < 0.2 ~ "E2+",
#                                     TRUE ~"Other")) |>
#         ggplot(aes(x = logFC, y = -log10(P.Value), color = DE_class)) +
#         geom_point(alpha = 0.5, size = 0.5) +
#         geom_text_repel(aes(label = ifelse(DE_class != "Other", Symbol, NA)), size = 2, show.legend=FALSE) +
#         scale_color_manual(values = c(, "Other" = "darkgray")) +
#         theme_bw() +
#         geom_vline(xintercept = c(1,0,-1), linetype = c("dashed", "solid","dashed")) +
#         geom_hline(yintercept = -log10(pval_lim), linetype = "dashed") +
#         labs(title = title, subtitle = subtitle)
# 
# 
# ggsave(custom_volcano, filename = here(plot_dir, "custom_volcano_carrier_model.png"), width = 12, height =4)

#### plot genes by carrier ####
plot_DE_express_carrier <- function(gene, cluster){
    
    df <- extractData(res.proc, assay = cluster, genes = gene)
    
    # expression boxplot
    expression_plot <- ggplot(df, aes(APOE_carrier, !!sym(gene), fill = APOE_carrier)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1) +
        ylab(bquote(Expression ~ (log[2] ~ CPM))) +
        labs(title = sprintf("%s - %s", cluster, gene)) +
        scale_fill_manual(values = APOE_carrier_colors) +
        theme_bw()
    
    ggsave(expression_plot, filename = here(plot_dir, sprintf("sn_carrier_expression_boxplot-%s-%s.png", cluster, gene)))
    
}

## boxplot all signif carrier DEGs
tt_signif |> 
    filter(model == "carrier") |>
    pmap(~plot_DE_express(cluster = ..1, gene = ..2))


#### forest plot ####
my_plot_forest <- function(gene, mod_name){
    forest_test <- plotForest(res.dl.list[[mod_name]], coef = coef_list[[mod_name]], gene = gene) + labs(subtitle = sprintf("mod=%s, coef=%s", mod_name, coef_list[[mod_name]]))
    ggsave(forest_test, filename = here(plot_dir, sprintf("dreamlet_sn_forest-%s-%s.png", mod_name, gene)))
}

## forest plot all signif carrier DEGs
tt_signif |> 
    filter(model == "carrier") |>
    pmap(~my_plot_forest(mod_name = "carrier", gene = ..2))

#### contrast DEGs ####
res.dl.contrast <- readRDS(here("processed-data", "10_dreamlet_sn", "04_run_dreamlet_contrast_sn", "dreamlet_contrast_sn-contrast.RDS"))
coefNames(res.dl.contrast)

contrast_coef <- c("E2E2_E4E4", 
                   "E3E4_E4E4", 
                   "E2E3_E4E4",
                   "E2E2_E3E4",
                   "E2E2_E2E3",
                   "E2E3_E3E4",
                   "E4E4_anyE2",
                   "anyE4_anyE2",
                   "E2E2_anyE4",
                   "E2E3_anyE4")

names(contrast_coef) <- contrast_coef

## get contrast top tables
tt_contrast <- map(contrast_coef, 
                   ~topTable(res.dl.contrast, coef = .x, number = Inf) 
)

## adjust pval
tt_contrast <- map2(tt_contrast, names(tt_contrast), 
                    ~.x  |>
                            group_by(assay) |>
                            mutate(adj.P.Val.cell_type = p.adjust(P.Value),
                                   contrast = .y)
                    )

## save
save(tt_contrast, file = here(data_dir, "dreamlet_sn_contrast_TopTables.Rdata"))

## summarize n signif
fdr_contrast_count <- map2_dfr(tt_contrast, names(tt_contrast), ~.x |> 
                                ungroup() |>
                                summarise(global_FDR10 = sum(adj.P.Val < 0.1),
                                          global_FDR20 = sum(adj.P.Val < 0.2),
                                          cell_type_FDR10 = sum(adj.P.Val.cell_type < 0.1),
                                          cell_type_FDR20 = sum(adj.P.Val.cell_type < 0.2)
                                          ) |>
                                mutate(model = .y) |>
                                as.data.frame()) 

write.csv(fdr_contrast_count, file = here(data_dir, "dreamlet_sn_FDR_contrast_count.csv"), row.names = FALSE)

## filter and export key results
tt_contrast_signif <- map_dfr(tt_contrast, ~.x |> as.data.frame() |> filter(adj.P.Val.cell_type < 0.2))

write.csv(tt_contrast_signif, file = here(data_dir, "dreamlet_sn_topTable_contrast_FDR20.csv"), row.names = FALSE)

#### Contrast volcano plot ####
pdf(here(plot_dir, "sn_dreamlet_VolcanoPlot_contrast.pdf"), height = 11, width = 8)

map(contrast_coef, ~plotVolcano(res.dl.contrast, .x, cutoff = 0.2) + labs(title = .x))

dev.off()

#### summarize contrast DEGs ####
tt_contrast_signif |> 
    group_by(assay, ID) |>
    summarise(contrasts = paste0(contrast, collapse = ", "))


tt_contrast_signif_summary <- tt_contrast_signif |> 
    group_by(assay, ID) |>
    summarise(contrasts = paste0(contrast, collapse = ", ")) |>
    arrange(ID) |>
    mutate(model = "contrast")
    
allDE_summary <- tt_signif |>
    group_by(assay, ID) |>
    summarise(model = paste0(model, collapse = ", ")) |>
    bind_rows(tt_contrast_signif_summary)

## Overlap w/ visium gene
visium_de_summary <- read.csv(here("processed-data", "11_dreamlet_Visium", "05_compile_dreamlet_Visium", "dreamlet_Visium_modelcontrast_summary.csv")) |>
    mutate(visium_DE = paste0(assay, ": ", models)) |>
    select(visium_DE, ID)

intersect(visium_de_summary$ID, allDE_summary$ID)

# visium_de_summary |> filter(ID == "CNTNAP4")

allDE_summary2 <- allDE_summary |>
    group_by(assay, ID) |>
    summarise(n_models = n(),
              models = paste0(unique(model), collapse = ", "),
              contrasts = paste0(contrasts[!is.na(contrasts)], collapse = ", ")) |>
    mutate(risk = ID %in% AD_risk$symbol) |>
    left_join(visium_de_summary)

allDE_summary2 |> filter(risk)

allDE_summary2 |> print(n= 49)

write.csv(allDE_summary2, file = here(data_dir, "dreamlet_sn_modelcontrast_summary.csv"), row.names = FALSE)


tt_contrast_signif |>
    dplyr::count(assay)
# assay  n
# 1 Astro 13
# 2 Excit 10
# 3 Inhib  5
# 4 Macro  9
# 5 Micro  1
# 6   OPC  6
# 7 Oligo  6
# 8  Vasc  7

#### plot expression by apoe ####
plot_DE_express_apoe <- function(gene, cluster, subtitle = NULL){
    
    df <- extractData(res.proc, assay = cluster, genes = gene)
    
    # expression boxplot
    expression_plot <- ggplot(df, aes(APOE, !!sym(gene), fill = APOE)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1) +
        ylab(bquote(Expression ~ (log[2] ~ CPM))) +
        labs(title = sprintf("%s - %s", cluster, gene), subtitle = subtitle) +
        scale_fill_manual(values = APOE_genotype_colors) +
        theme_bw() 
    
    ggsave(expression_plot, filename = here(plot_dir, sprintf("sn_apoe_expression_boxplot-%s-%s.png", cluster, gene)))
    
}

pmap(tt_contrast_signif_summary, function(...) plot_DE_express_apoe(cluster = ..1, gene = ..2, subtitle = ..3))

#### plot expression by apoe + ancestry ####
plot_DE_express_apoe_anc <- function(gene, cluster, subtitle = NULL){
    
    df <- extractData(res.proc, assay = cluster, genes = gene)
    
    # expression boxplot
    expression_plot <- ggplot(df, aes(APOE, !!sym(gene), fill = APOE, color = Ancestry)) +
        geom_boxplot() +
        # geom_jitter(width = .1) +
        ylab(bquote(Expression ~ (log[2] ~ CPM))) +
        labs(title = sprintf("%s - %s", cluster, gene), subtitle = subtitle) +
        scale_fill_manual(values = APOE_genotype_colors) +
        scale_color_manual(values = ancestry_colors) +
        theme_bw() 
    
    ggsave(expression_plot, filename = here(plot_dir, sprintf("sn_apoe_anc_expression_boxplot-%s-%s.png", cluster, gene)))
    
}

pmap(tt_signif |> filter(model == "apoe_i"), function(...) plot_DE_express_apoe_anc(cluster = ..1, gene = ..2, subtitle = "apoe_i"))

#### contrast vs. regular model ####
compare_contrast_scatter <- function(mod, contrast){
    
    contrast_compare <- tt[[mod]] |>
        as.data.frame() |>
        select(assay, ID, mod_logFC = logFC, mod_t = t, mod_pval = adj.P.Val.cell_type) |>
        full_join(tt_contrast[[contrast]] |>
                      as.data.frame() |>
                      select(assay, ID, contrast_logFC = logFC, contrast_t = t, contrast_pval = adj.P.Val.cell_type)) |>
        mutate(DE_class = case_when(mod_pval < 0.2 & contrast_pval < 0.2 ~ "signif_both",
                                    mod_pval < 0.2 ~ "signif_mod",
                                    contrast_pval < 0.2 ~ "signif_contrast",
                                    TRUE ~"None"
        ))
    
    # contrast_compare |> count(DE_class)
    # summary(contrast_compare)
    
    contrast_compare_scatter <- contrast_compare |>
        ggplot(aes(x = mod_t, y = contrast_t, color = DE_class)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_text_repel(aes(label = ifelse(DE_class != "None", ID, "")), size = 1.5) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = c(signif_both = "purple", signif_mod = "blue", signif_contrast = "red")) +
        theme_bw() +
        facet_wrap(~assay) +
        labs(title = sprintf("%s vs. %s", mod, contrast))
    
    ggsave(contrast_compare_scatter, filename = here(plot_dir, sprintf("sn_contrast_compare_scatter_%s-v-%s.png", mod, contrast)), height = 10, width = 10)
    
}

compare_contrast_scatter("carrier", "anyE4_anyE2")
compare_contrast_scatter("carrier", "E4E4_anyE2")
compare_contrast_scatter("e4e4", "E2E2_E4E4")
compare_contrast_scatter("e4e4", "E4E4_anyE2")

#### Bernie DE ####

# readRDS("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/12_DEanalyses_removedsampsAndFinalNMseg/02c-Clusterwise_DETopTabs_Sex_APO_predomAncest.RDS")


# slurmjobs::job_single('05_compile_dreamlet_sn', create_shell = TRUE, memory = '10G', command = "Rscript 05_compile_dreamlet_sn.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


