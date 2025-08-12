## Louise Huuki-Myers, Aug 2025
## Compile LDSC results
## based on https://github.com/LieberInstitute/spatialdACC/blob/main/code/17-2_LDSC_enrichment/spatial/ldsc_results2.R

#### Set Up ####

library("tidyverse")
library("getopt")
library("here")
library("sessioninfo")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type",
      "mode", "m", "1", "character", "specificity or enrichment"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)


opt$datatype <- "sn_broad"
opt$mode <- "specificity"

data_dir <- here("processed-data", "18_LDSC", "04_ldsc_results")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Get output ####
out_dir <- here("processed-data", "18_LDSC", sprintf("LDSC_%s_%s", opt$datatype, opt$mode))
clusters <- list.files(out_dir)

ldsc_results <- map_dfr(clusters, function(clus){
    
    ldsc_resuls_files <- list.files(here(out_dir, clus), pattern = ".results", full.names = TRUE)
    gwas <- gsub(".out.results", "", basename(ldsc_resuls_files))
    ldsc_results <- map_dfr(ldsc_resuls_files, ~read.table(.x, header = TRUE) |> filter(Category == "L2_0")) |>
        add_column(cluster = clus, gwas = gwas,  .before = 1) 
    
    return(ldsc_results)
    
})

head(ldsc_results)

#### adjust p-vale ####

ldsc_results <- ldsc_results |>
    mutate(p_zcore =  pnorm(abs(Coefficient_z.score),lower.tail=F)*2,
           FDR = p.adjust(p_zcore, method="fdr"))

## save data

write.csv(ldsc_results, file = here(data_dir, sprintf("LDSC_results-%s_%s.scv", opt$datatype, opt$mode)))

# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "01_prep_marker_bed",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



