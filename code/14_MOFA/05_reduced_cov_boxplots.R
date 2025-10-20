library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
plot_path = here("plots", "14_MOFA", "05_reduced_cov_boxplots", "boxplot.pdf")
specific_factor = "Factor3"

dir.create(dirname(plot_path), showWarnings = FALSE)

model = load_model(model_path)
