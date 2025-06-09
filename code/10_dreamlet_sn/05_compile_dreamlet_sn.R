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
data_dir <- here("processed-data", "10_dreamlet_sn", "04_plot_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "04_plot_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))
# 
# res.dl <- readRDS(here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn", sprintf("dreamlet_sn-%s.RDS", opt$model)))

res.dl.fn <- list.files(here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn"), full.names = TRUE)
names(res.dl.fn) <- gsub(".RDS", "", gsub("dreamlet_sn-", "", basename(res.dl.fn)))

res.dl.list <- purrr::map(res.dl.fn, readRDS)
purrr::map(res.dl.list, coefNames)

coef_list <- list(apoe_i = c("APOEE2/E3:Anc_Afr", "APOEE3/E4:Anc_Afr", "APOEE4/E4:Anc_Afr"),
                  apoe_n0 = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  apoe = c("APOEE2/E3", "APOEE3/E4","APOEE4/E4"),
                  carrier_i = "APOE_carrierE4+:Anc_Afr",
                  carrier_n0 = "APOE_carrierE4+",
                  carrier_n1 = "APOE_carrierE4+",
                  carrier_n2 = "APOE_carrierE4+",
                  carrier_n3 = "APOE_carrierE4+",
                  carrier_n4 = "APOE_carrierE4+",
                  carrier_sf = "APOE_carrierE4+",
                  carrier = "APOE_carrierE4+",
                  contrast = c("APOE_synE2.E2","APOE_synE2.E3","APOE_synE3.E4","APOE_synE4.E4"),
                  e4e4_i = "APOE_E4E4TRUE:Anc_Afr",
                  e4e4_n0 = "APOE_E4E4TRUE",
                  e4e4 = "APOE_E4E4TRUE")

identical(names(coef_list), names(res.dl.list))

details(res.dl.list[["carrier"]])

topTable(res.dl.list[["carrier"]], "APOE_carrierE4+")

topTable(res.dl.list[["contrast"]], c("APOE_synE2.E2", "APOE_synE4.E4"))

# results from full analysis
tt <- map2(res.dl.list, coef_list, 
           ~topTable(.x, .y, number = Inf, p.value = 0.1) |> 
               group_by(assay) |>
               mutate(adj.P.Val.cell_type = p.adjust(P.Value))
)

fdr_01 <- map(tt, ~.x |> as.data.frame() |> filter(adj.P.Val < 0.1))
map_int(fdr_01, nrows)

map(tt, ~.x |> as.data.frame() |> group_by(assay) |> filter(adj.P.Val < 0.2) |> count())

map(tt, ~.x |> summarise(global_signigf = sum(adj.P.Val < 0.2),
                         cell_type_signif = sum(adj.P.Val.cell_type < 0.2)))

map2_dfr(tt, names(tt), ~.x |> 
        ungroup() |>
        summarise(global_signigf = sum(adj.P.Val < 0.1),
                  cell_type_signif = sum(adj.P.Val.cell_type < 0.1)) |>
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



pdf(here(plot_dir, "sn_dreamlet_VolcanoPlot.pdf"), height = 11, width = 8)
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

ggsave(expression_plot, filename = here(plot_dir, "sn_dreamlet_expression_boxplot.png"))

## forest plot
plotForest(res.dl, coef = "APOE_carrierE4+", gene = "NBPF12")

# slurmjobs::job_single('04_plot_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 04_plot_dreamlet_sn.R")

# slurmjobs::job_loop(
#     loops = list(model = names(dreamlet_models_sn)),
#     name = "04_plot_dreamlet_sn",
#     create_shell = TRUE,
#     create_script = FALSE
# )

