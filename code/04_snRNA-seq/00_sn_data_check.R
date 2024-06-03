
library("tidyverse")
library("readxl")
library("here")

metadata_sn <- read.csv(here("processed-data", "00_project_prep", "02_get_online_metadata","metadata_sn_plan.csv"))
metadata_sn

sample_info_fn <- list.files(here("raw-data", "sample_info_snRNAseq"), full.names = TRUE)
sample_info <- map_dfr(sample_info_fn, read_excel)

summary(sample_info)

# Nuclei Sorted   Nuclei Targeted cDNA Amp Cycle Agilent [cDNA] pg/ul Dilution Factor...11 Final [cDNA] pg/ul
# Min.   :14000   Min.   :9000    Min.   :11     Min.   : 435.4       Min.   :1.000        Min.   : 435.4    
# 1st Qu.:14000   1st Qu.:9000    1st Qu.:11     1st Qu.:1158.1       1st Qu.:1.000        1st Qu.:1412.5    
# Median :14000   Median :9000    Median :11     Median :1472.8       Median :1.000        Median :1812.1    
# Mean   :14000   Mean   :9000    Mean   :11     Mean   :1635.9       Mean   :1.214        Mean   :1863.3    
# 3rd Qu.:14000   3rd Qu.:9000    3rd Qu.:11     3rd Qu.:2003.8       3rd Qu.:1.000        3rd Qu.:2263.9    
# Max.   :14000   Max.   :9000    Max.   :11     Max.   :3304.7       Max.   :2.000        Max.   :3304.7    
# 
# Total cDNA ng      cDNA Input       SI cycles     Ave frag length Agilent [lib] pg/ul Dilution Factor...18
# Min.   : 17.42   Min.   : 4.354   Min.   :14.00   Min.   :438.0   Min.   :1091        Min.   : 4.000      
# 1st Qu.: 56.50   1st Qu.:14.125   1st Qu.:14.00   1st Qu.:453.2   1st Qu.:2000        1st Qu.: 6.000      
# Median : 72.48   Median :18.121   Median :14.00   Median :460.5   Median :2143        Median : 8.500      
# Mean   : 74.53   Mean   :18.633   Mean   :14.57   Mean   :470.4   Mean   :2096        Mean   : 9.964      
# 3rd Qu.: 90.56   3rd Qu.:22.639   3rd Qu.:15.00   3rd Qu.:486.5   3rd Qu.:2431        3rd Qu.:12.250      
# Max.   :132.19   Max.   :33.047   Max.   :16.00   Max.   :523.0   Max.   :2569        Max.   :24.000      
# 
# Final [lib] pg/ul  index_name         index(i7)         index2_workflow_a(i5) index2_workflow_b(i5)
# Min.   : 2355     Length:14          Length:14          Length:14             Length:14            
# 1st Qu.: 9228     Class :character   Class :character   Class :character      Class :character     
# Median :14837     Mode  :character   Mode  :character   Mode  :character      Mode  :character     
# Mean   :17912                                                                                      
# 3rd Qu.:21473                                                                                      
# Max.   :47348                                                                                      
# 
# Est Read Pairs (million)
# Min.   :420             
# 1st Qu.:420             
# Median :420             
# Mean   :420             
# 3rd Qu.:420             
# Max.   :420 

# cn <- map(sample_info, colnames)
# table(unlist(cn))

(cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger")))
# [1] "10c_ERC_SVB" "11c_ERC_SVB" "12c_ERC_SVB" "13c_ERC_SVB" "14c_ERC_SVB" "1c_ERC_SVB"  "4c_ERC_SVB"  "5c_ERC_SVB" 
# [9] "6c_ERC_SVB"  "7c_ERC_SVB"  "8c_ERC_SVB"  "9c_ERC_SVB" 

sample_info |>
    select(sample_id = `Sample #`, BrNum = Brain, exp_round = `Experimental Round`, seq_round = `Sequencing Round`) |>
    mutate(path = here("processed-data", "03_cellranger", sample_id, "outs", "raw_feature_bc_matrix"),
           file_exists = file.exists(path)
           )
# sample_id   BrNum  Round Seq_Round done 
# <chr>       <chr>  <dbl>     <dbl> <lgl>
# 1 1c_ERC_SVB  Br1039     1         1 TRUE 
# 2 2c_ERC_SVB  Br1556     1        NA FALSE
# 3 3c_ERC_SVB  Br5212     1        NA FALSE
# 4 4c_ERC_SVB  Br2582     2         1 TRUE 
# 5 5c_ERC_SVB  Br5276     2         1 TRUE 
# 6 6c_ERC_SVB  Br5415     2         1 TRUE 
# 7 7c_ERC_SVB  Br1691     3         1 TRUE 
# 8 8c_ERC_SVB  Br3974     3         1 TRUE 
# 9 9c_ERC_SVB  Br5426     3         1 TRUE 
# 10 10c_ERC_SVB Br6476     3         1 TRUE 
# 11 11c_ERC_SVB Br1706     4         1 TRUE 
# 12 12c_ERC_SVB Br2305     4         1 TRUE 
# 13 13c_ERC_SVB Br5460     4         1 TRUE 
# 14 14c_ERC_SVB Br5517     4         1 TRUE 

slurmjobs::job_loop(
    loops = list(sample_id = cell_ranger_out_paths),
    name = "01_get_droplet_scores",
    create_shell = TRUE
)

