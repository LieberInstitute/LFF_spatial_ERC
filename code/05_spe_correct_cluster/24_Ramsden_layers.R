## Louise Huuki-Myers, April 2025
## Compare expression to DEGs from Ramsden et al., 2015 

library("spatialLIBD")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("Orthology.eg.db")

## set up dirs
plot_dir <- here("plots", "05_spe_correct_cluster", "24_Ramsden_layers")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "24_Ramsden_layers")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load Ramsden Gene lists ####
# Ensembl mouse database (release 73)

ramsden_sheets <- c("L2only", "L3only", "L5only")
names(ramsden_sheets) <- ramsden_sheets
ramsden_data <- map(ramsden_sheets, ~read_excel(here("external-data", "Ramsden2015", "pcbi.1004032.s011.xls"), sheet = .x))

## other annotations
ramsden_sheets2 <- c("LayerPatternedLayer5a6all",'LayerPatternedLayer6all', "LayerPatternedLayer5ball","LayerPatternedIslandsPara","LayerPatternedIslands","LayerPatternedL2b")
names(ramsden_sheets2) <- ramsden_sheets2
ramsden_data_other <- map(ramsden_sheets2, ~read_excel(here("external-data", "Ramsden2015", "pcbi.1004032.s011.xls"), sheet = .x))


## All Gene sheet
ramsden_data_all <- read_excel(here("external-data", "Ramsden2015", "pcbi.1004032.s011.xls"), sheet = "AllGene")

ramsden_data_all |> count(Vis) |> print(n = 32)
# filter to layer specific - drop None, Uniform, error ect. 
ramsden_data_all <- ramsden_data_all |> filter(grepl("^L", Vis))

ramsden_data_list <- map(splitit(ramsden_data_all$Vis), ~ramsden_data_all$Name[.x])
map_int(ramsden_data_list, length)

## all the same as individual lists? - Yes
map2_lgl(ramsden_data, names(ramsden_data), ~setequal(.x$Name, ramsden_data_all |> filter(Vis == .y) |> pull(Name)))

## look up emsembl w/ entrez
mouse_entrez_to_human <- function(mouse_entrez){
    mouse_entrez <- mouse_entrez[mouse_entrez != "0"]
    mapped <- select(Orthology.eg.db, mouse_entrez, "Homo_sapiens","Mus_musculus")
    hu_ensembl <- select(org.Hs.eg.db, as.character(mapped[,2]), "ENSEMBL","ENTREZID")
    return(hu_ensembl)
}
#test
# mouse_entrez_to_human(c("66643", "16498"))

ramsden_data_list_entrez <- map(splitit(ramsden_data_all$Vis), ~ramsden_data_all$Entrez[.x])
ramsden_data_list_entrez <- c(ramsden_data_list_entrez, map(ramsden_data_other, "Entrez"))
ramsden_data_human <- map(ramsden_data_list_entrez, ~mouse_entrez_to_human(.x))

ramsden_human_genes <- map(ramsden_data_human, "ENSEMBL")
map_int(ramsden_human_genes, length)

# L232    L233 L23even    L252    L255 L25even  L2high   L2mid  L2only  L2weak    L355 L35even  L3high   L3mid 
# 45      22     290      42      29     169     371     709      76     124       5      24      59      31 
# L3only  L5high   L5mid  L5only  L5weak 
# 31      83     241      56      19 

#### gene enrichment ####

## get reference layer enrichment statistics
# 
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

## enrichment

ramsden_enrichment <- gene_set_enrichment(
    gene_list = ramsden_human_genes,
    modeling_results = erc_modeling,
    model_type = "enrichment"
)

#### plot enrichment ####
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

ramsden_enrichment <- ramsden_enrichment |> mutate(test = factor(test, levels = names(SpD_colors)))

pdf(here(plot_dir, "ramsden_ernichment.pdf"), width = 12, height = 8)
gene_set_enrichment_plot(enrichment = ramsden_enrichment,
                         plot_SetSize_bar = TRUE,
                         model_colors = SpD_colors
)
dev.off()

ramsden_enrichment_only <- ramsden_enrichment |> filter(grepl("only", ID))

pdf(here(plot_dir, "ramsden_ernichment_only.pdf"), width = 12, height = 8)
gene_set_enrichment_plot(enrichment = ramsden_enrichment_only,
                         plot_SetSize_bar = TRUE,
                         model_colors = SpD_colors
)
dev.off()
