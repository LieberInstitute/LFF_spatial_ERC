## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")
library("HDF5Array")

#### Set up dirs ####
data_dir <- here("processed-data", "11_dreamlet_Visium", "01_prep_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "11_dreamlet_Visium", "01_prep_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

## Use gene names
rownames(spe) <- rowData(spe)$gene_name
## deal with duplicated gene_names
rownames(spe)[duplicated(rownames(spe))] <- rowData(spe)$gene_id[duplicated(rownames(spe))]
stopifnot(!any(duplicated(rownames(spe))))

## add syntactic APOE vars
spe$APOE_syn <- factor(gsub("/", ".", spe$APOE))
levels(spe$APOE_syn)

spe$APOE_carrier_syn <- factor(gsub("\\+", "", spe$APOE_carrier))
levels(spe$APOE_carrier_syn)

## add E4/E4 variable
spe$APOE_E4E4 <- spe$APOE == "E4/E4"
table(spe$APOE_E4E4)

## syntatic SpD
spe$SpD_syn <- gsub("~", "_", spe$SpD)

## drop uneeded 10x clusters -> lead to error
colData(spe)[,grep("^10x", colnames(colData(spe)))] <- NULL

#### Aggregate to pseudobulk ####
message(Sys.time(), ' - Aggregate to pseudobulk')
pb <- aggregateToPseudoBulk(spe,
                            assay = "counts",
                            cluster_id = "SpD_syn",
                            sample_id = "sample_id",
                            verbose = FALSE
)

## one assay per SpD
assayNames(pb)

dim(pb)

# Normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, ~APOE, min.count = 5)

# the resulting object of class dreamletProcessedData stores
# normalized data and other information
res.proc

details(res.proc)

# view details of dropping samples
details(res.proc)

## save 
saveRDS(res.proc, file = here(data_dir, "Visium_res_proc.rds"))

# res.proc <- readRDS(here("processed-data", "11_dreamlet_Visium", "01_prep_dreamlet_sn", "Visium_res_proc.rds"))

# show voom plot for each cell clusters
pdf(here(plot_dir, "Visium_dreamlet_voom.pdf"))
plotVoom(res.proc)
dev.off()


metadata(res.proc)

form <- ~ APOE + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium)

model.matrix(~0+APOE_syn, colData(res.proc))

# L <- makeContrastsDream(~0+APOE, colData(res.proc),
#                         contrasts = c(
#                             `E2/E2-E4/E4` = "APOEE2/E2 - APOEE4/E4",
#                             `E3/E4-E4/E4` = "APOEE3/E4 - APOEE4/E4",
#                             `E2/E3-E4/E4` = "APOEE2/E3 - APOEE4/E4",
#                             `E2/E2-E3/E4` = "APOEE2/E2 - APOEE3/E4",
#                             `E2/E2-E2/E3` = "APOEE2/E2 - APOEE2/E3",
#                             `E2/E3-E3/E4` = "APOEE2/E3 - APOEE3/E4")
# )

L <- makeContrastsDream(~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium_slide), colData(res.proc),
                        contrasts = c(
                            `E2E2_E4E4` = "APOE_synE2.E2 - APOE_synE4.E4",
                            `E3E4_E4E4` = "APOE_synE3.E4 - APOE_synE4.E4",
                            `E2E3_E4E4` = "APOE_synE2.E3 - APOE_synE4.E4",
                            `E2E2_E3E4` = "APOE_synE2.E2 - APOE_synE3.E4",
                            `E2E2_E2E3` = "APOE_synE2.E2 - APOE_synE2.E3",
                            `E2E3_E3E4` = "APOE_synE2.E3 - APOE_synE3.E4")
)

# Visualize contrast matrix
contrast_plot <- plotContrasts(L) + labs(title = "DREAM contrast",
                                         subtitle = "~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium_slide)")

ggsave(contrast_plot, filename = here(plot_dir, "Visium_dream_contrast.png"))

# slurmjobs::job_single('01_prep_dreamlet_Visium', create_shell = TRUE, memory = '25G', command = "Rscript 01_prep_dreamlet_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

