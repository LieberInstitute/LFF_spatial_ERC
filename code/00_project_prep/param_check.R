
library("getopt")

spec <- matrix(
    c(  "cluster", "i", "1", "character", "Name of cluster",
        "k", "k", "2", "numeric", "Number of clusters",
        "drop_WM", "d", "3", "logical", "drop White Matter - defined in script"
    ),
    ncol = 5, byrow = TRUE
)

opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

opt <- list(cluster ="BayesSpace_SVGm", k = 10, drop_WM = TRUE)

lapply(opt, "class")

# Rscript param_check.R --cluster "BayesSpace_SVGm" --k $10 --drop_WM TRUE