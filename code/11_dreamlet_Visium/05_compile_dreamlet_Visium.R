## Louise Huuki-Myers, June 2025
## Compile and plot dreamlet sn data

library("SpatialExperiment")
library("dreamlet")
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("DFplyr")

# Import command-line parameters
# scec <- matrix(
#     c("model", "m", "1", "character", "Model name"),
#     ncol = 5, byrow = TRUE
# )
# opt <- getopt(scec)
# print(opt)

# test 
# opt$model <- "carrier"

#### Set up dirs ####
data_dir <- here("processed-data", "11_dreamlet_Visium", "05_compile_dreamlet_Visium")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "11_dreamlet_Visium", "05_compile_dreamlet_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- readRDS(here("processed-data", "11_dreamlet_Visium", "01_prep_dreamlet_Visium", "Visium_res_proc.rds"))
# 
# res.dl <- readRDS(here("processed-data", "11_dreamlet_Visium", "03_run_dreamlet_Visium", sprintf("dreamlet_Visium-%s.RDS", opt$model)))

res.dl.fn <- list.files(here("processed-data", "11_dreamlet_Visium", "03_run_dreamlet_Visium"), full.names = TRUE)
names(res.dl.fn) <- gsub(".RDS", "", gsub("dreamlet_Visium-", "", basename(res.dl.fn)))

res.dl.list <- purrr::map(res.dl.fn, readRDS)
res.dl.list[["contrast"]] <- NULL
purrr::map(res.dl.list, coefNames)

coef_list <- list(apoe_i = c("APOEE2/E3:Anc_Afr", "APOEE3/E4:Anc_Afr", "APOEE4/E4:Anc_Afr"),
                  apoe_n0 = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  apoe = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  carrier_i = "APOE_carrierE4+:Anc_Afr",
                  carrier_kr = "APOE_carrierE4+",
                  carrier_n0 = "APOE_carrierE4+",
                  carrier_n1 = "APOE_carrierE4+",
                  carrier_n2 = "APOE_carrierE4+",
                  carrier_n3 = "APOE_carrierE4+",
                  carrier_n4 = "APOE_carrierE4+",
                  carrier_sf = "APOE_carrierE4+",
                  carrier = "APOE_carrierE4+",
                  # contrast = c("APOE_synE2.E2","APOE_synE2.E3","APOE_synE3.E4","APOE_synE4.E4"),
                  e4e4_i = "APOE_E4E4TRUE:Anc_Afr",
                  e4e4_n0 = "APOE_E4E4TRUE",
                  e4e4 = "APOE_E4E4TRUE")

identical(names(coef_list), names(res.dl.list))

details(res.dl.list[["carrier"]])

topTable(res.dl.list[["carrier"]], "APOE_carrierE4+")

## sex DGE
topTable(res.dl.list[["carrier_sf"]], "SexM")

topTable(res.dl.list[["carrier_sf"]], "SexM", number = Inf) |>
    group_by(assay) |>
    mutate(adj.P.Val.cell_type = p.adjust(P.Value)) |>
    summarise(global_signigf = sum(adj.P.Val < 0.1),
              cell_type_signif = sum(adj.P.Val.cell_type < 0.1))

# groups[xx, -ncol(groups)] global_signigf cell_type_signif
# <character>      <integer>        <integer>
# 1             Inhib_Sp09D09              6                2
# 2                L1_Sp09D05              4                3
# 3              L2.3_Sp09D01              6                2
# 4                L5_Sp09D03              7                5
# 5                L6_Sp09D04             10                8
# 6                LD_Sp09D02              4                2
# 7              Vasc_Sp09D08              6                5
# 8                WM_Sp09D06              7                3
# 9             WM.uf_Sp09D07              0                0


# results from full analysis
tt <- map2(res.dl.list, coef_list, 
           ~topTable(.x, .y, number = Inf, p.value = 0.1) |> 
               group_by(assay) |>
               mutate(adj.P.Val.SpD = p.adjust(P.Value))
)

fdr10 <- map(tt, ~.x |> filter(adj.P.Val < 0.1))

map(tt, ~.x |> filter(adj.P.Val.SpD < 0.1))
map_int(fdr10, nrow)

map(tt, ~.x |> as.data.frame() |> group_by(assay) |> filter(adj.P.Val < 0.2) |> count())

map(tt, ~.x |> summarise(global_signigf = sum(adj.P.Val < 0.1),
                         SpD_signif = sum(adj.P.Val.SpD < 0.1)))

map2_dfr(tt, names(tt), ~.x |> 
        ungroup() |>
        summarise(global_signigf = sum(adj.P.Val < 0.1),
                  SpD_signif = sum(adj.P.Val.SpD < 0.1)) |>
         mutate(model = .y) |>
            as.data.frame())



## check VarPart genes
tt[["carrier"]] |>
    as.data.frame() |>
    filter(ID %in% c("CLOCK", "ABCC4"))

tt[["apoe"]] |>
    as.data.frame() |>
    filter(ID %in% c("ADAM10", "ADAM17"))

tt[["e4e4"]] |>
    as.data.frame() |>
    filter(ID %in% c("SRGAP2C", "NAALADL2", "APP"))



pdf(here(plot_dir, "Visium_dreamlet_VolcanoPlot.pdf"), height = 11, width = 8)
# plotVolcano(res.dl, coef = "APOE_carrierE4+")

map(c("carrier","e4e4","carrier_i","e4e4_i"), ~plotVolcano(res.dl.list[[.x]], coef_list[[.x]]) + labs(title = .x))

dev.off()

#### plot genes ####
# get data
df <- extractData(res.proc, "Vasc", genes = "TAFA1")

# expression boxplot
expression_plot <- ggplot(df, aes(APOE_carrier, TAFA1)) +
    geom_boxplot() +
    ylab(bquote(Expression ~ (log[2] ~ CPM))) +
    ggtitle("NBPF12") +
    theme_bw()

ggsave(expression_plot, filename = here(plot_dir, "Visium_dreamlet_expression_boxplot.png"))

## forest plot
plotForest(res.dl, coef = "APOE_carrierE4+", gene = "NBPF12")

# slurmjobs::job_single('05_compile_dreamlet_Visium', create_shell = TRUE, memory = '50G', command = "Rscript 05_compile_dreamlet_Visium.R")

# slurmjobs::job_loop(
#     loops = list(model = names(dreamlet_models_Visium)),
#     name = "05_compile_dreamlet_Visium",
#     create_shell = TRUE,
#     create_script = FALSE
# )

