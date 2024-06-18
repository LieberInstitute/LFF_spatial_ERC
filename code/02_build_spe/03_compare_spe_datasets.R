
library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "02_build_spe", "02_explore_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### LIBD spatial data ####

## access with spatialLIBD

spatialLIBD_datasets <- c(visium_pilot = "spe", 
                          spatialDLPFC = "spatialDLPFC_Visium",
                          visium_AD = "Visium_SPG_AD_Visium_wholegenome_spe")

spatial_data_libd <- map2_dfr(spatialLIBD_datasets, names(spatialLIBD_datasets), 
                         function(lookup, name){
                            spe <- spatialLIBD::fetch_data(lookup)
                            cn <- colnames(colData(spe))
                            print(cn)
                            qc_cols <- c("sample_id", "in_tissue", "sum_umi", "sum_gene", "expr_chrM_ratio")
                            qc <- as.data.frame(colData(spe)[,qc_cols[qc_cols %in% cn]])
                            qc$dataset <- name
                            return(qc)
                           }
                         )

map(spatial_data_libd, head)
spatial_data_libd |>
    count(dataset)

## access on JHPCE
jhpce_datasets <- c(
    LFF_spatial_ERC = here("processed-data", "02_build_spe", "spe.rds"),
    LFF_spatial_ERC_trim = here("processed-data", "02.1_spe_compare_trim", "spe_trim.rds"),
    spatial_NAc = "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/05_harmony_BayesSpace/01-build_spe/spe_raw.rds",
    # spatial_HYP = "/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/processed-data/02_build_spe/spe.Rdata",
    LFF_spatial_LC = "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/02_build_spe/spe.rds",
    Habenula_Visium = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/02_build_spe/spe.rds"
    )

map(jhpce_datasets, file.exists)

spatial_data_jhpce <- map2_dfr(jhpce_datasets, names(jhpce_datasets), 
                              function(path, name){
                                  message(Sys.time(), " - loading ", name)
                                  spe <- readRDS(path)
                                  cn <- colnames(colData(spe))
                                  # print(cn)
                                  qc_cols <- c("sample_id", "in_tissue", "sum_umi", "sum_gene", "expr_chrM_ratio")
                                  qc <- as.data.frame(colData(spe)[,qc_cols[qc_cols %in% cn]])
                                  qc$dataset <- name
                                  return(qc)
                              }
)

spatial_data_qc <- bind_rows(spatial_data_jhpce, spatial_data_libd)
spatial_data_qc |> count(in_tissue, dataset)

sample_qc <- spatial_data_qc |>
    group_by(dataset, sample_id) |>
    filter(in_tissue) |>
    summarize(sum_umi = median(sum_umi),
              sum_gene = median(sum_gene),
              expr_chrM_ratio = median(expr_chrM_ratio),
              )|>
    pivot_longer(!c(dataset, sample_id), names_to = "metric", values_to = "median")

median_metic_boxplots <- sample_qc |>
    ggplot(aes(x = dataset, y = median)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    facet_wrap(~metric, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(median_metic_boxplots, filename = here(plot_dir, "median_metric_boxplots.png"), width = 10)


