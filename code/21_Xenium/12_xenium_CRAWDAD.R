## Louise Huuki-Myers, May 2026
## Run CRAWDAD on xenium data

library(crawdad)
library(tidyverse)
## load the spleen data of the pkhl sample 
data('pkhl')
## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], cellTypes = pkhl$celltypes)
## define the scales to analyze the data
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
## shuffle cells to create null background
shuffle_list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = 7,
                                            seed = 1,
                                            verbose = TRUE)
## calculate the zscore for the cell-type pairs at different scales
results <- crawdad::findTrends(cells,
                               neighDist = 50,
                               shuffleList = shuffle_list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization
vizColocDotplot(dat, zSigThresh = zsig, zScoreLimit = 2*zsig, 
                dotSizes = c(3,15)) +
    theme(axis.text.x = element_text(angle = 35, h = 0))