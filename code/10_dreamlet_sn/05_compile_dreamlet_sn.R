## Louise Huuki-Myers, June 2025
## Compile and plot dreamlet sn data

library("SpatialExperiment")
library("dreamlet")
library("tidyverse")
library("here")
library("sessioninfo")
library("DFplyr")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "05_compile_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "05_compile_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

#### Load data ####
res.proc <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))
# 
# res.dl <- readRDS(here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn", sprintf("dreamlet_sn-%s.RDS", opt$model)))

res.dl.fn <- list.files(here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn"), full.names = TRUE)
names(res.dl.fn) <- gsub(".RDS", "", gsub("dreamlet_sn-", "", basename(res.dl.fn)))

res.dl.list <- purrr::map(res.dl.fn, readRDS)
purrr::map(res.dl.list, coefNames)

coef_list <- list(apoe_i = c("APOEE2/E3:Anc_Afr", "APOEE3/E4:Anc_Afr", "APOEE4/E4:Anc_Afr"),
                  apoe_n0 = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  apoe = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  carrier_i = "APOE_carrierE4+:Anc_Afr",
                  carrier_n0 = "APOE_carrierE4+",
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
tt <- map2(res.dl.list, coef_list, 
           ~topTable(.x, .y, number = Inf) |> 
               group_by(assay) |>
               mutate(adj.P.Val.cell_type = p.adjust(P.Value))
)

tt <- map2(tt, names(tt), ~.x |> mutate(model = .y))


map(tt, ~.x |> summarise(global_FDR10 = sum(adj.P.Val < 0.1),
                         cell_type_FDR10 = sum(adj.P.Val.cell_type < 0.1)))

fdr_model_count <- map2_dfr(tt, names(tt), ~.x |> 
        ungroup() |>
        summarise(global_FDR10 = sum(adj.P.Val < 0.1),
                  cell_type_FDR10 = sum(adj.P.Val.cell_type < 0.1)) |>
         mutate(model = .y) |>
            as.data.frame()) 

write.csv(fdr_model_count, file = here(data_dir, "dreamlet_sn_FDR10_model_count.csv"), row.names = FALSE)

## filter and export key results
main_mods <- c("carrier", "apoe", "e4e4", "carrier_i", "apoe_i", "e4e4_i")

tt_signif <- map_dfr(tt[main_mods], ~.x |> as.data.frame() |> filter(adj.P.Val.cell_type < 0.2))

write.csv(tt_signif, file = here(data_dir, "dreamlet_sn_topTable_FDR20.csv"), row.names = FALSE)

save(tt[main_mods], file = here(data_dir, "dreamlet_sn_TopTables.Rdata"))

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
        geom_histogram(binwidth = 0.1) +
        facet_wrap(~assay) +
        theme_bw() +
        labs(title = mod)
    
    ggsave(pval_histo, filename = here(plot_dir, sprintf("dreamlet_Vsium_pval_histo-%s.png", mod)), width = 10)
    
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

#### Volvano Plots ####
pdf(here(plot_dir, "sn_dreamlet_VolcanoPlot.pdf"), height = 11, width = 8)
# plotVolcano(res.dl, coef = "APOE_carrierE4+")

map(c("carrier","e4e4","carrier_i","e4e4_i"), ~plotVolcano(res.dl.list[[.x]], coef_list[[.x]]) + labs(title = .x))

dev.off()

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

tt_contrast <- map(contrast_coef, 
                   ~topTable(res.dl.contrast, coef = .x, number = Inf) 
)

## save
save(tt_contrast, file = here(data_dir, "dreamlet_sn_contrast_TopTables.Rdata"))

tt_contrast <- map2(tt_contrast, names(tt_contrast), 
                    ~.x  |>
                            group_by(assay) |>
                            mutate(adj.P.Val.cell_type = p.adjust(P.Value),
                                   contrast = .y)
                    )

fdr_contrast_count <- map2_dfr(tt_contrast, names(tt_contrast), ~.x |> 
                                ungroup() |>
                                summarise(global_FDR10 = sum(adj.P.Val < 0.1),
                                          cell_type_FDR10 = sum(adj.P.Val.cell_type < 0.1)) |>
                                mutate(model = .y) |>
                                as.data.frame()) 

write.csv(fdr_contrast_count, file = here(data_dir, "dreamlet_sn_FDR10_contrast_count.csv"), row.names = FALSE)

## filter and export key results
tt_contrast_signif <- map_dfr(tt_contrast, ~.x |> as.data.frame() |> filter(adj.P.Val.cell_type < 0.1))

write.csv(tt_contrast_signif, file = here(data_dir, "dreamlet_sn_topTable_contrast_FDR10.csv"), row.names = FALSE)

#### Contrast volcano plot ####
pdf(here(plot_dir, "sn_dreamlet_VolcanoPlot_contrast.pdf"), height = 11, width = 8)

map(contrast_coef, ~plotVolcano(res.dl.contrast, .x) + labs(title = .x))

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

visium_de_summary <- read.csv(here("processed-data", "11_dreamlet_Visium", "05_compile_dreamlet_Visium", "dreamlet_Visium_modelcontrast_summary.csv")) 
    # select(visium_de = assay, ID)


allDE_summary2 <- allDE_summary |>
    group_by(assay, ID) |>
    summarise(n_models = n(),
              models = paste0(unique(model), collapse = ", "),
              contrasts = paste0(contrasts[!is.na(contrasts)], collapse = ", ")) |>
    mutate(risk = ID %in% AD_risk$symbol) 

allDE_summary2 |> filter(risk)

allDE_summary2 |> print(n= 29)

## no overlap w/ visium genes :(
intersect(visium_de_summary$ID, allDE_summary2$ID)

write.csv(allDE_summary2, file = here(data_dir, "dreamlet_sn_modelcontrast_summary.csv"), row.names = FALSE)


tt_contrast_signif |>
    dplyr::count(assay)

#   assay n
# 1 Astro 5
# 2 Excit 4
# 3 Inhib 2
# 4 Macro 2
# 5   OPC 3
# 6 Oligo 2
# 7  Vasc 4

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


# slurmjobs::job_single('05_compile_dreamlet_sn', create_shell = TRUE, memory = '10G', command = "Rscript 05_compile_dreamlet_sn.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


