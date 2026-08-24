
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
                      # Macro = "#79354E",
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

APOE_carrier_colors_dark <- c(`E2+`="#3B918B", `E4+`="#B7522A")
APOE_carrier_colors <- c(`E2+`="#51B8B1", `E4+`="#D97D59")

## phenotype colors
ancestry_colors <- c(EA="#8BA1E4",AA="#A4C77F")
sex_colors <- c(M = "#5C80BC", F ="#D58BCC")

## APOE ancestry colors 

APOE_carrier_anc_colors <- c(
    "E2_EA" = "#51B8B1", 
    "E2_AA" = "#A8D9D6", 
    "E4_EA" = "#D97D59",    
    "E4_AA" = "#EDB99E"  
)

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
save(ancestry_colors, 
     sex_colors, 
     APOE_genotype_colors, 
     APOE_carrier_colors, 
     APOE_carrier_colors_dark,
     sample_colors,
     APOE_carrier_anc_colors,
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

#### Cell Type Colors V2 ####
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

cell_type_colors_broad.V2 <- c(Astro = "#3BB273",
                               Macro = "#79354E",
                               Micro = "#663894",
                               Oligo = "#F57A00",
                               OPC = "#D2B037",
                               Vasc = "#FF56AF",
                               Excit = "#247FBC",
                               Inhib = "#E83E38",
                               Other = "#4E586A")

## plot official colors
pdf(here(plot_dir, "ERC_cell_type_colors_broad.V2.pdf"))
test_pallet_plots(cell_type_colors, "ERC Broad Cell Type Colors V2")
dev.off()

cell_type_colors_anno = c(Astro.1 = "#228B22",
                          Astro.2 = "#32CD32",
                          Astro.3 = "#808000",
                          Astro.4 = "#98FF98",
                          Astro.5 = "#50C878",
                          Macro = "#79354E",
                          Micro.1 = "#663894",
                          Micro.2 = "#A375D1",
                          Micro.3 = "#B96FD9",
                          Micro.4 = "#673AB7",
                          Micro.5 = "#BFA8ED",
                          OPC.1 = "#D2B037",
                          OPC.2 = "#BDB76B",
                          OPC.3 = "#FFDB58",
                          OPC.4 = "#A2852D",
                          OPC.5 = "#DA9100",
                          Oligo.1 = "#F57A00",
                          Oligo.2 = "#F2A240",
                          Oligo.3 = "#CC5500",
                          Oligo.4 = "#C27153",
                          Oligo.5 = "#FF8C00",
                          Vasc.Endo = "#FF56AF",
                          Vasc.PC = "#DE289E",
                          Vasc.VLMC = "#DA5FE8", 
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

save(cell_type_colors_anno, file = here("processed-data","00_project_prep","cell_type_colors_anno_subtype.Rdata"))

cell_type_colors <- list(broad = cell_type_colors_broad.V2,
                         anno = cell_type_colors_anno)

save(cell_type_colors, file = here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"))

## Oligo + OPC alt colors
# Oligo_OPC_colors <- create_cell_colors(cell_types = c(paste0("OPC.", 1:5), paste0("Oligo.", 1:5)), palette_name = "gg")
# OPC.1     OPC.2     OPC.3     OPC.4     OPC.5   Oligo.1   Oligo.2   Oligo.3   Oligo.4   Oligo.5
# "#F8766D" "#D89000" "#A3A500" "#39B600" "#00BF7D" "#00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC"

Oligo_OPC_colors <- c(cell_type_colors_anno[grepl("OPC", names(cell_type_colors_anno))],
                      Oligo.1 = "#00BFC4",
                      Oligo.2 =  "#00B0F6",
                      Oligo.3 = "#9590FF",
                      Oligo.4 = "#E76BF3",
                      Oligo.5 = "#FF62BC")

# OPC.1     OPC.2     OPC.3     OPC.4     OPC.5   Oligo.1   Oligo.2   Oligo.3   Oligo.4   Oligo.5 
# "#D2B037" "#BDB76B" "#FFDB58" "#A2852D" "#DA9100" "#46aed7" "#6169d6" "#976bb3" "#7088cd" "#b25fce" 

save(Oligo_OPC_colors, file = here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"))

#### SpD colors ####
# SpD_colors_V1 <- c("Vasc~Sp09D08" = "#E05AD2", #Orchid
#                 "L1~Sp09D05" = "#0220DE", #Chrystler Blue
#                 "L2.3~Sp09D01" = "#FEAF16", #light orange
#                 "LD~Sp09D02" = "#00BCF9", #dark sky blue
#                 "Inhib~Sp09D09" = "#C82100", #Engineering red
#                 "L5~Sp09D03" = "#16FF32", #lime
#                 "L6~Sp09D04" = "#178C6D", #forest green
#                 "WM.uf~Sp09D07" = "#E4E1E3", # purple white
#                 "WM~Sp09D06" = "#581009") #brown

SpD_colors <- c(
    "vVasc"      = "#E05AD2",
    "vL1"        = "#16C72B",
    "vL2"        = "#021AB6",
    "vInhib"     = "#C82100",
    "vL3"        = "#889DF0",
    "vLD"        = "grey70",
    "vL5"        = "#0087F5",
    "vL6"        = "#40DAF2",
    "vWMuf"      = "#F4A460",
    "vWMim"      = "#E8720C",
    "vWMd"       = "#581009"
)


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


SpD_colors_V3 <- c("Vasc" = "#E05AD2", 
                "L1a" = "#16FF32", 
                "L1b" = "#178C6D",
                "L2.3" = "#FEAF16", 
                "LD" = "#5A5156",
                "Inhib" = "#C82100", 
                "L5" = "#00BCF9", 
                "L6" = "#0220DE", 
                "WM.uf" = "#E4E1E3",
                "WM" = "#581009") 


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

#### Spatial DLPFC k9 colors ####
spatialDLPFC_SpD_colors <- c(`L1~Sp09D01` = "#5A5156",
                             `L1~Sp09D02` =   "#E4E1E3",
                             `L2~Sp09D03` =  "#F6222E",
                             `L3~Sp09D05` =  "#16FF32",
                             `L4~Sp09D08` = "#B00068",
                             `L5~Sp09D04` = "#FE00FA",
                             `L6~Sp09D07` =  "#FEAF16",
                             `WM~Sp09D06`  = "#3283FE", 
                             `WM~Sp09D09` = "#1CFFCE")
         
save(spatialDLPFC_SpD_colors, file = here("processed-data","00_project_prep", "spatialDLPFC_Data","spatialDLPFC_SpD_colors.Rdata"))




## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()