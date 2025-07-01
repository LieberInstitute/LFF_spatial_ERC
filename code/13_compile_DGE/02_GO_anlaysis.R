## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"
# opt$datatype = "Visium"

data_dir <- here("processed-data", "13_compile_DGE", "02_GO_analysis")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "02_GO_analysis")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))

## load DE data ##
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn)

DE_data |> count(cluster)

head(DE_data)
DE_data |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster)

## ENTREZID look up
entrez_search <- bitr(DE_data$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

DE_entrez <- DE_data |> 
    left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "up",
                                vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class)) ## doesn't like .  in cluster names
    
DE_entrez |> count(cluster)
DE_entrez |> filter(DE_class != "None") |> count(DE_class_cluster)

# DE_clusters <- DE_entrez |> filter(DE_class != "None") |> pull(cluster) |> unique()
# names(DE_clusters) <- DE_clusters

universe <- unique(DE_entrez$ENTREZID)
length(universe)

ont <- c("CC","BP","MF")
names(ont) <- ont

go_result <- map(ont, ~compareCluster(ENTREZID ~ DE_class_cluster,
                                      data = DE_entrez |> filter(DE_class != "None"), 
                                      OrgDb = org.Hs.eg.db,
                                      fun = enrichGO,
                                      universe = universe,
                                      ont = .x, ##ALL,CC,BP,MF
                                      pAdjustMethod = "BH",
                                      # pvalueCutoff = 1,#0.5,
                                      qvalueCutoff = 0.05,
                                      readable = TRUE))


pdf(file = here(plot_dir, sprintf("GO_dotplot_%s.pdf", opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), 
      ~print(
          dotplot(.x, 
                  x = "DE_class_cluster", 
                  showCategory = 5, 
                  label_format = 60)  +
              ggtitle(paste("GO Enrichment:", .y)) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
      )
      )
dev.off()

compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y))

compare_clus |> count(DE_class_cluster)


write.csv(compare_clus, file = here(data_dir, sprintf("GO_results_%s.csv", opt$datatype)), row.names = FALSE)

# slurmjobs::job_single('02_GO_analysis', create_shell = TRUE, memory = '5G', command = "Rscript 02_GO_analysis.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

