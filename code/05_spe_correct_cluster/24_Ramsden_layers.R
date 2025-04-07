## Louise Huuki-Myers, April 2025
## Compare expression to DEGs from Ramsden et al., 2015 

library("spatialLIBD")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")

# library("cowplot")
# library("patchwork")

library("org.Hs.eg.db")
# library("RCurl")
library("AnnotationHub")

## source reduced dims function

plot_dir <- here("plots", "05_spe_correct_cluster", "24_Ramsden_layers")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "24_Ramsden_layers")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load Ramsden Gene lists ####
# Ensembl mouse database (release 73)

# ramsden_sheets <- c("L2only", "L3only", "L5only")
# names(ramsden_sheets) <- ramsden_sheets
# ramsden_data <- map(ramsden_sheets, ~read_excel(here("external-data", "Ramsden2015", "pcbi.1004032.s011.xls"), sheet = .x))

ramsden_data_all <- read_excel(here("external-data", "Ramsden2015", "pcbi.1004032.s011.xls"), sheet = "AllGene")

ramsden_data_all |> count(Vis) |> print(n = 32)

ramsden_data_all <- ramsden_data |> filter(grepl("^L", Vis))

ramsden_data_list <- map(splitit(ramsden_data_all$Vis), ~ramsden_data_all$Name[.x])

# Basic function to convert mouse to human gene names from https://support.bioconductor.org/p/9153600/
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

mouse_to_human <- function(mouseids, horg, morg, orth){
    mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
    mapped <- select(orth, mouseg, "Homo_sapiens","Mus_musculus")
    names(mapped) <- c("Mus_egid","Homo_egid")
    husymb <- select(horg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
    # huensembl <- select(horg, as.character(mapped[,2]), "ENSEMBL","ENTREZID")
    return(data.frame(Mus_symbol = mouseids,
                      mapped,
                      Homo_symbol = husymb[,2]))
}

# ramsden_data_human <- map(ramsden_data, ~mouse_to_human(.x$Name, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db))
ramsden_data_human <- map(ramsden_data_list, ~mouse_to_human(.x, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db))

# TODO debug ENSEMBL lookup
# mouse_to_human_ensembl <- function(mouseids, horg, morg, orth){
#     mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
#     mapped <- select(orth, mouseg, "Homo_sapiens","Mus_musculus")
#     names(mapped) <- c("Mus_egid","Homo_egid")
#     husymb <- select(horg, as.character(mapped[,2]), "ENSEMBL","ENTREZID")
#     return(data.frame(Mus_symbol = mouseids,
#                       mapped,
#                       Homo_symbol = husymb[,2]))
# }
# 
# mouse_to_human_ensembl(c('Tle1', "Grp"), org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)

# ramsden_data_human <- map(ramsden_data, ~mouse_to_human_ensembl(.x$Name, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db))
# ℹ With name: L2only.
# Caused by error in `data.frame()`:
#     ! arguments imply differing number of rows: 76, 63

ramsden_human_genes <- map(ramsden_data_human, ~unique(.x$Homo_symbol[!is.na(.x$Homo_symbol)]))
map_int(ramsden_human_genes, length)
# L2only L3only L5only 
# 58     21     46

#### gene enrichment ####

## get reference layer enrichment statistics

# erc_modeling <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm", "modeling_results-BayesSpace_SVGm_k09.rds"))
# 
# erc_spd_anno <- readxl::read_excel(here("processed-data","05_spe_correct_cluster", "10_spatial_registration_DLPFC", "ERC_SpD_spatial_registration_Annotations.xlsx")) |>
#     mutate(Annotation = fct_reorder(Annotation, order)) |>
#     mutate(SpD = fct_reorder(paste0(Annotation, "~", cluster), order)) |>
#     dplyr::select(cluster, Annotation, SpD)
# 
# erc_modeling <- map(erc_modeling, function(mod){
#     erc_colnames <- colnames(mod)
#     pwalk(erc_spd_anno, function(...) erc_colnames <<- gsub(..1, ..3, erc_colnames))
#     colnames(mod) <- erc_colnames
#     return(mod)
# })

erc_modeling <- readRds(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "modeling_results-BayesSpace_SVGm_k09_annotated.rds"))

erc_modeling$enrichment$gene

table(duplicated(erc_modeling$enrichment$gene))

map_lgl(ramsden_human_genes, ~"HSPA14" %in% .x)

erc_modeling$enrichment <- erc_modeling$enrichment[erc_modeling$enrichment$gene != "HSPA14",]

rownames(erc_modeling$enrichment) <- erc_modeling$enrichment$gene
erc_modeling$enrichment$ensembl <- erc_modeling$enrichment$gene


## enrichment

ramsden_enrichment <- gene_set_enrichment(
    gene_list = ramsden_human_genes,
    modeling_results = erc_modeling,
    model_type = "enrichment"
)



pdf(here(plot_dir, "ramsden_ernichment.pdf"))
gene_set_enrichment_plot(enrichment = ramsden_enrichment,
                         plot_SetSize_bar = TRUE)
dev.off()
