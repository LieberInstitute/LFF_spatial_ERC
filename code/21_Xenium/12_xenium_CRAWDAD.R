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

# spe <- spe[, spe$spot_class == "singlet"]'

spe <- spe[, spe$BrNum == opt$brnum]

#### covert data ####

#   Gather spatial coordinates for input to CRAWDAD
pod_tb <- tibble(x = spatialCoords(spe)[, 'x_centroid'],
                 y = spatialCoords(spe)[, 'y_centroid'])

## convert to sf
cells <- crawdad:::toSF(pos = spatial_cells_sample[, c("x", "y")],
                        cellTypes = spe$first_type)


# pdf(here(plot_dir, "viz_cluster_test.pdf"))
# crawdad::vizEachCluster(cells = cells,
#                         coms = as.factor(spe$first_type),
#                         s = 2)
# dev.off()

## define the scales to analyze the data
scales = c(100, 200, 500, 1000, 5000)
random_seed = 513

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))


## shuffle cells to create null background
shuffle_list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = num_cores,
                                            seed = random_seed,
                                            verbose = TRUE)

## calculate the zscore for the cell-type pairs at different scales
results <- crawdad::findTrends(cells,
                               neighDist = 50,
                               shuffleList = shuffle_list,
                               ncores = num_cores,
                               verbose = TRUE,
                               returnMeans = FALSE)

dat <- crawdad::meltResultsList(results, withPerms = TRUE)

## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization

vizColocDotplot(dat, zSigThresh = zsig, zScoreLimit = 2*zsig, 
                dotSizes = c(3,15)) +
    theme(axis.text.x = element_text(angle = 35, h = 0))


