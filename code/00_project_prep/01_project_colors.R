
library("colorblindr")
library("tidyverse")
library("spatialLIBD")
library("here")

source(here("code", "utils", "test_pallet_plots.R"))

plot_dir <- here("plots", "00_project_prep", "01_project_colors")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### cell types ####
cell_type_colors <- c(Excit = "#247FBC",
                      Inhib = "#E83E38",
                      Oligo = "#F57A00",
                      OPC = "#D2B037",
                      Astro = "#3BB273",
                      Micro = "#663894",
                      Macro = "#79354E",
                      Endo = "#FF56AF",
                      # PC = "#DE289E",
                      # VLMC = "#DA5FE8",
                      Other = "#4E586A")

## plot official colors
pdf(here(plot_dir, "ERC_cell_type_colors.pdf"))
test_pallet_plots(cell_type_colors, "Cell Colors: ERC")
dev.off()

save(cell_type_colors, file = here("processed-data","00_project_prep","cell_type_colors.Rdata"))

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)

## APOE colors
# APOE_genotype_colors <- c(`E2/E2`="#114B5F", `E2/E3` = "#61C9A8", `E3/E4`="#ED9B40", `E4/E4`="#BA3B46")
APOE_genotype_colors <- c(`E2/E2`="#186E8B", `E2/E3` = "#61C9A8", `E3/E4`="#ED9B40", `E4/E4`="#BA3B46")
APOE_carrier_colors <- c(`E2+`="#398A84", `E4+`="#D46B43")

## phenotype colors
ancestry_colors <- c(EA="#1B3174",AA="#698F3F")
sex_colors <- c(M = "#5C80BC", F ="#D58BCC")

## test APOE colors
pdf(here(plot_dir, "ERC_APOE_colors.pdf"))
test_pallet_plots(APOE_genotype_colors, "APOE_genotype_colors")
dev.off()

#### sample colors ####
sample_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

apoe_colors_tb <- sample_info |> 
    select(sample_id, APOE) |> 
    unique() |>
    group_by(APOE) |>
    arrange(sample_id) |>
    mutate(APOE_color = paste0(APOE, ".", row_number())) |>
    arrange(APOE_color)

apoe_colors <- DeconvoBuddies::create_cell_colors(apoe_colors_tb$APOE_color, pallet = APOE_genotype_colors, split = "\\.")
sample_colors <- apoe_colors$fine
names(sample_colors) <- apoe_colors_tb$sample_id

# Br1556    Br2582    Br5529    Br5712    Br5832    Br6161    Br5212    Br5367    Br5415    Br5426    Br5634    Br5854    Br6423
# "#186E8B" "#3E869E" "#659EB1" "#8BB6C4" "#B2CED8" "#D8E6EB" "#61C9A8" "#74CFB2" "#88D6BD" "#9CDDC8" "#B0E3D3" "#C3EADE" "#D7F1E9"
# Br6538    Br3974    Br6476    Br5276    Br5460    Br5517    Br5599    Br6085    Br6098    Br6263    Br6321    Br1039    Br1289
# "#EBF8F4" "#ED9B40" "#EEA453" "#F0AF66" "#F2B979" "#F4C38C" "#F6CD9F" "#F7D7B2" "#F9E1C5" "#FBEBD8" "#FDF5EB" "#BA3B46" "#C35760"
# Br1691    Br1706    Br2305    Br5161    Br5941
# "#CD727A" "#D78F95" "#E1ABAF" "#EBC6CA" "#F5E3E4"

## save
save(ancestry_colors, sex_colors, APOE_genotype_colors, APOE_carrier_colors, sample_colors,
     file = here("processed-data", "project_colors.Rdata"))

#### Cell type color archive ####

## other cell type colors to test
#
# ## remix for ERC
# cell_type_remix <- c(Excit = "#374194",
#                       Inhib = "#f50018", 
#                       Oligo = "#f57a00", 
#                       OPC = "#d1cf36",
#                       Astro = "#3BB273",
#                       Micro = "#7a3894",
#                       Endo = "#e8389b",
#                       Other = "#4E586A")
# 
# cell_type_colors_alt <- c(Excit = "#246EBD",
#                       Inhib = "#E83870",
#                       Oligo = "#DD900B",
#                       OPC = "#BAD136",
#                       Astro = "#3BB391",
#                       Micro = "#4C3894",
#                       Endo = "#FA63C8",
#                       Other = "#4E586A")
# 
# cell_type_colors_cool <- c(Excit = "#2423BC",
#                            Inhib = "#94378A",
#                            Oligo = "#E8A738",
#                            OPC = "#DDF500",
#                            Astro = "#8CC82D",
#                            Micro = "#3BA9B2",
#                            Endo = "#5C537F")
# 
# cell_type_colors_list <- list(classic = cell_type_colors,
#                               remix = cell_type_remix,
#                               alt = cell_type_colors_alt,
#                               cool = cell_type_colors_cool)

# walk2(cell_type_colors_list, names(cell_type_colors_list), function(colors, name){
#     pdf(here(plot_dir, sprintf("cell_type_color_%s_test_plots.pdf", name)))
#     test_pallet_plots(colors, paste0("Cell Colors: ", name))
#     dev.off()
# })

#### more cell colors ####
list(Crimson = "#DC143C",
     Firebrick = "#B22222",
     Indian_Red = "#CD5C5C",
     Light_Coral = "#F08080",
     Dark_Red = "#8B0000",
     Tomato = "#FF6347")

list(Red_Violet = "#C71585",
     Magenta_Red = "#D21868",
     Crimson_Red = "#DC143C",
     Pure_Red = "#FF0000",
     Scarlet_Red = "#FF4500",
     Red_Orange = "#FF6347")

list(Royal_Blue = "#4169E1",
     Dodger_Blue = "#1E90FF",
     Steel_Blue = "#4682B4",
     Sky_Blue = "#87CEEB",
     Midnight_Blue = "#191970",
     Cornflower_Blue = "#6495ED")

list(Green_Blue = "#00B3B3",
     Teal_Blue = "#0099CC",
     Cerulean_Blue = "#007ACC",
     Azure_Blue = "#005FCC",
     Cobalt_Blue = "#0040CC",
     Violet_Blue = "#2E2EFF")

cell_type_colors_anno = c(Astro.1 = "#3BB273",
                          Astro.2 = "#9DD8B9",
                          Endo = "#FF56AF",
                          Micro.1 = "#663894",
                          Micro.2 = "#B29BC9",
                          Oligo.1 = "#F57A00",
                          Oligo.2 = "#F8A655",
                          choird_plexis = "#AA7F6E",
                          OPC = "#D2B037",
                          Excit.L2="#4169E1",
                          Excit.L2_5.1="#1E90FF",
                          Excit.L2_5.2="#4682B4",
                          Excit.L5.1 = "#87CEEB",
                          Excit.L5.2="#212197",
                          Excit.L5_6_NP="#6495ED",
                          Excit.L6_CT="#0099CC",
                          Excit.L6b="#2E2EFF",
                          Inhib.Pax6="#DC143C",
                          Inhib.Lamp5_Lhx6="#B22222",
                          Inhib.Pvalb="#CD5C5C",
                          Inhib.Vip="#F08080",
                          Inhib.Chandelier="#E83E38",
                          Inhib.Sst="#8B0000")

pdf(here(plot_dir, "ERC_cell_type_colors_anno.pdf"), height = 15, width = 11)
test_pallet_plots(cell_type_colors_anno, "Cell Colors Anno: ERC")
dev.off()

save(cell_type_colors_anno, file = here("processed-data","00_project_prep","cell_type_colors_anno.Rdata"))


#### SpD colors ####
SpD_colors <- c("Vasc~Sp09D08" = "#E05AD2", #Orchid
                "L1~Sp09D05" = "#0220DE", #Chrystler Blue
                "L2.3~Sp09D01" = "#FEAF16", #light orange
                "L3~Sp09D02" = "#00BCF9", #dark sky blue
                "LD~Sp09D09" = "#C82100", #Engineering red
                "L5~Sp09D03" = "#16FF32", #lime
                "L6~Sp09D04" = "#178C6D", #forest green
                "WM.uf~Sp09D07" = "#E4E1E3", # purpe white
                "WM~Sp09D06" = "#581009") #brown

pdf(here(plot_dir, "ERC_SpD_colors.pdf"), height = 11, width = 8)
test_pallet_plots(SpD_colors, "SpD Colors: ERC")
dev.off()

save(SpD_colors, file = here("processed-data","SpD_colors.Rdata"))

SpD_colors <- c("Vasc~Sp09D08" = "#E05AD2",
                "L1~Sp09D05" = "#0220DE",
                "L2.3~Sp09D01" = "#80428A",
                "L3~Sp09D02" = "#AFADFF",
                "LD~Sp09D09" = "#F55200",
                "L5~Sp09D03" = "#A2E838",
                "L6~Sp09D04" = "#147B5F",
                "WM.uf~Sp09D07" = "#E4E1E3",
                "WM~Sp09D06" = "#581009")

save(SpD_colors, file = here("processed-data","00_project_prep","SpD_colors2.Rdata"))

## plot offical colors
pdf(here(plot_dir, "ERC_SpD_colors.2.pdf"), height = 11, width = 8)
test_pallet_plots(SpD_colors, "SpD Colors: ERC2")
dev.off()

## SpD archive colors
# "#16FF32"    "#90AD1C"    "#5A5156"    "#3283FE"      "brown"    "#FE00FA"    "#F6222E"    "#FEAF16"    "#1CFFCE"

## archive
# c("Spanish orange"="#E96F00",
# "jonquil yellow"="#FFCD17",
# "Citrine"="#E3D348",
# "Phthalo blue"="#021380",
# "Black bean"="#500802")



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()