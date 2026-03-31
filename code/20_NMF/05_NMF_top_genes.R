## Louise Huuki-Myers, Mar 2026
## Select top genes for factors

#### Set up ####

# library("singlet")
# library("dplyr")
# library("ggplot2")

library("SingleCellExperiment")
library("RcppML")
library("UpSetR")
library("scDotPlot")
library("sessioninfo")
library("here")
library(tidyverse)

data_dir <- here("processed-data", "20_NMF", "05_NMF_top_genes")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "20_NMF", "05_NMF_top_genes")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### load data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

## Drop Br1289
sce <- sce[, sce$BrNum != "Br1289"]

rd <- rowData(sce) |>
    as.data.frame() |>
    select(gene_name, gene_id) |>
    as_tibble()

## NMF data
nmf <- readRDS(here("processed-data", "20_NMF", "01_runNMF", "ERC_sn_nmf.RDS"))

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(nmf, decreasing=TRUE)[1:n]))
    return(top_genes)
}

apply(nmf@w, 2, function(x) names(sort(x, decreasing=TRUE)[1:10]))


factor_weights_long <- nmf@w |>
    reshape2::melt() |>
    rename(gene_id = Var1, Factor = Var2, weight = value) |>
    group_by(Factor) |>
    mutate(weight_rank = rank(-weight)) |>
    left_join(rd)

factor_top_genes <- factor_weights_long |>
    arrange(weight_rank) |>
    slice(1:10)

write.csv(factor_top_genes, file = here(data_dir, "ERC_sn_nmf_top_genes.csv"))

## overlap?

factor_overlap <- factor_top_genes |>
    group_by(gene_name) |>
    summarise(n = n(),
           factors = list(Factor))

factor_overlap |> arrange(-n)
#    gene_name     n factors   
#     <chr>     <int> <list>    
#     1 MALAT1       17 <fct [17]>
#     2 NRXN1        11 <fct [11]>
#     3 CNTNAP2       9 <fct [9]> 
#     4 CSMD1         9 <fct [9]> 
#     5 RBFOX1        9 <fct [9]> 
#     6 DPP10         7 <fct [7]> 
#     7 ERBB4         7 <fct [7]> 
#     8 KCNIP4        7 <fct [7]> 
#     9 PCDH9         7 <fct [7]> 
#     10 IL1RAPL1      6 <fct [6]> 

factor_overlap |> count(n) |> arrange(-n)

## create named list 
factor_top_genes_list <- map(rafalib::splitit(factor_top_genes$Factor), ~factor_top_genes$gene_name[.x])

ncol(nmf@h) == ncol(sce)


