## Louise Huuki-Myers, Aug 2025
## Find overlap in DEG and LR data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("data.table")

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

data_dir <- here("processed-data", "17_cell_cell_comm", "10_LR_DEG_check")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "17_cell_cell_comm", "10_LR_DEG_check")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# load(here("processed-data", "project_colors.Rdata"))

#### load DE data ####
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn) |> select(gene_name, cluster, vlmf_t, vlmf_adj.P.Val, vlmf_logFC)

DE_data_signif <- DE_data |> filter(vlmf_adj.P.Val < 0.05)
DE_data_signif |> count(cluster)



#top pairs 
LR_top_pairs <- tibble(LR = c('NPTX1^NPTXR', 'RTN4^LINGO1', 'GNAS^ADCY1', 'RTN4^RTN4R', 'CALM1^PTPRA', 'PSAP^SORT1')) |>
    separate(LR, sep= "\\^", into = c("L", "R"), remove = FALSE)

LR_top_pairs_long <- LR_top_pairs |> 
    pivot_longer(!LR, names_to = "class", values_to = "gene_name")
    
LR_top_pairs_long |> left_join(DE_data_signif)
# LR          class gene_name cluster vlmf_t vlmf_adj.P.Val vlmf_logFC
# <chr>       <chr> <chr>     <fct>    <dbl>          <dbl>      <dbl>
# 1 NPTX1^NPTXR L     NPTX1     NA       NA           NA          NA    
# 2 NPTX1^NPTXR R     NPTXR     Astro.1   5.73         0.0310      0.751
# 3 NPTX1^NPTXR R     NPTXR     Oligo.3   4.16         0.0178      1.21 
# 4 RTN4^LINGO1 L     RTN4      NA       NA           NA          NA    
# 5 RTN4^LINGO1 R     LINGO1    NA       NA           NA          NA    
# 6 GNAS^ADCY1  L     GNAS      NA       NA           NA          NA    
# 7 GNAS^ADCY1  R     ADCY1     NA       NA           NA          NA    
# 8 RTN4^RTN4R  L     RTN4      NA       NA           NA          NA    
# 9 RTN4^RTN4R  R     RTN4R     Oligo.3   3.18         0.0443      1.48 
# 10 CALM1^PTPRA L     CALM1     NA       NA           NA          NA    
# 11 CALM1^PTPRA R     PTPRA     NA       NA           NA          NA    
# 12 PSAP^SORT1  L     PSAP      NA       NA           NA          NA    
# 13 PSAP^SORT1  R     SORT1     NA       NA           NA          NA   

#### bivariate data ####
bivariate_data <- fread(here("processed-data","17_cell_cell_comm","liana","bivariate_stats.csv.gz"))

bivariate_data |> filter(ligand == "ADGRB1", receptor == "RTN4R")

bivariate_data |> group_by(ligand, receptor) |> summarise()

# top-expressed and highest-spatial-association LR pairs

#### load LR data ####
liana_fn <- here("processed-data", "17_cell_cell_comm", "liana", "ranked_results", sprintf("%s.csv.gz", tolower(gsub("sn_", "", opt$datatype))))
file.exists(liana_fn)

liana_data <- fread(liana_fn)


