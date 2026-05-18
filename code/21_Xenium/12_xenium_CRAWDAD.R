## Louise Huuki-Myers, May 2026
## Run CRAWDAD on xenium data
## Adapted from https://github.com/LieberInstitute/Habenula_Visium/blob/b8719fbcdde7b14be267a396c1dd743cdc034649/code/09_HD_cell_level/15_crawdad_run.R

#### Set Up ####

library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("crawdad")
library("tidyverse")
library("getopt")
library("BiocParallel")

# Import command-line parameters
scec <- matrix(
    c("brnum", "b", "1", "character", "BrNum of selected sample"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# opt$brnum <- "Br1556"

data_dir <- here("processed-data", "21_Xenium", "12_xenium_CRAWDAD")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "12_xenium_CRAWDAD")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### load data ####
message(Sys.time(), " - Load SPE data, subset to: ", opt$brnum)
spe <- qs_read(here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_cell_types.qs2"))

spe <- spe[, spe$BrNum == opt$brnum]

message("ncells:", ncol(spe))

table(spe$cell_type_anno)
table(spe$cell_type_broad)

#### covert data ####

#   Gather spatial coordinates for input to CRAWDAD * NEEDS to be data.fame not tibble
coord_tb <- data.frame(x = spatialCoords(spe)[, 'x_centroid'],
                 y = spatialCoords(spe)[, 'y_centroid'])

## convert to sf
cells <- crawdad:::toSF(pos = coord_tb,
                        cellTypes = spe$cell_type_anno)

head(cells)

cells$celltypes <- droplevels(cells$celltypes)

levels(cells$celltypes)
table(cells$celltypes)

# pdf(here(plot_dir, "viz_cluster_test.pdf"))
# crawdad::vizEachCluster(cells = cells,
#                         coms = as.factor(spe$cell_type_broad),
#                         s = 2)
# dev.off()

## define the scales to analyze the data
scales = c(50, 100, 200, 300, 400, 500, 1000, 5000)
random_seed = 513

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# ## SET UP MULTICORE PARAM - not sure if needed
# ncores <- 4  # match Sys.getenv('SLURM_CPUS_ON_NODE')
bp <- MulticoreParam(workers = num_cores) #or bp <- MulticoreParam(4)


message("Run on ", num_cores, "cores")

## shuffle cells to create null background
shuffle_list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = num_cores,
                                            seed = random_seed,
                                            verbose = TRUE
                                            )

## calculate the zscore for the cell-type pairs at different scales
results <- crawdad::findTrends(cells, 
                               shuffleList = shuffle_list, 
                               returnMeans = FALSE,
                               ncores = num_cores, 
                               verbose = TRUE
)

dat <- crawdad::meltResultsList(results, withPerms = TRUE)

write_csv(dat |> mutate(BrNum = opt$brnum, .before = 1), file = here(data_dir, sprintf("xenium_CRAWDAD_results_%s.csv", opt$brnum)))

## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization

pdf(here(plot_dir, sprintf("xenium_CRAWDAD_colloc_dotplot_%s.pdf", opt$brnum)), height = 8, width = 10)
vizColocDotplot(dat, zSigThresh = zsig, zScoreLimit = 2*zsig, 
                dotSizes = c(3,15)) +
    theme(axis.text.x = element_text(angle = 35, h = 0))
dev.off()


#### slurm job ####
# xenium_experiment_info <- read.csv(here("processed-data", "21_Xenium", "07_xenium_build_spe", "xenium_experiment_details.csv"))
# 
# slurmjobs::job_loop(
#     loops = list(brnum = xenium_experiment_info$BrNum),
#     name = "12_xenium_CRAWDAD",
#     create_shell = TRUE,
#     create_script = FALSE
# )

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
