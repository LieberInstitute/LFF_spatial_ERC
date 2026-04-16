## Louise Huuki-Myers, April 2026
## prep xenium ranger run - NO NEED TO RUN matrix already in cell_feature_matrix/

library(here)

dir.create(here("processed-data", "21_Xenium", "07_xenium_ranger"))
# "/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/21_Xenium/07_xenium_ranger

xenium_files <- list.files(here("raw-data", "xenium"))
"/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/raw-data/xenium"



xenium_files <- data.frame(BrNum = gsub(".*(Br[0-9]{4}).*", "\\1", xenium_files),
                           dir = xenium_files)

write.table(xenium_files, file = here("code",  "21_Xenium", "xenium_samples.tsv"), col.names = FALSE, row.names = FALSE, quote = FALSE)


slurmjobs::job_single('07_xenium_ranger', create_shell = TRUE, memory = '25G', command = "xeniumranger")
