## Louise Huuki-Myers, April 2026
## Create summary plots all DE models

#### Set up ####
library("tidyverse")
library("here")

## set plot dir
data_dir <- here("plots", "13_compile_DGE", "16_all_mod_summary")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## set plot dir
plot_dir <- here("plots", "13_compile_DGE", "16_all_mod_summary")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### Get summary numbers ####

datatypes <- c("Visium", "sn_broad", "sn_fine")
names(datatypes) <- datatypes

mod_dir <- c(carrier = "01_compile_DGE",
             carrier_anc = "05_compile_DGE_ancestry")

# list.files(here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", paste0("vlmf_","sn_fine")), pattern = "FDR05_summary")

summary_files <- map2(mod_dir, names(mod_dir), function(mod, mod_name){
    unlist(map(c("Visium", "sn_broad", "sn_fine"), ~list.files(here("processed-data", "13_compile_DGE", mod, .x),
                                            pattern = "summary",
                                            full.names = TRUE))
           ) })

summary_interaction_tbl <- map2_dfr(summary_files, names(summary_files), 
                                    function(files, name){map(files, ~read.csv(.x, row.names = 1) |> 
                                                                  mutate(model = name,
                                                                         file = basename(.x),
                                                                         datatype = gsub(".*_summary_(.*)\\.csv", "\\1", file))
                                    )}
) |>
    filter(mod == "carrier" | is.na(mod))

summary_interaction_tbl |> count(model, datatype)

# summary_tbl <- map2_dfr(mod_dir, names(mod_dir), function(mod, mod_name){
#     map_dfr(datatypes, ~read.csv(list.files(here("processed-data", "13_compile_DGE", mod, paste0("vlmf_",.x)), 
#                                             pattern = "summary", 
#                                             full.names = TRUE), row.names = 1) |>
#                 mutate(Mod = mod_name, 
#                        datatype = .x))
# }) 


# summary_interaction_files <- map(datatypes, ~list.files(here("processed-data", "13_compile_DGE", "05_compile_DGE_interaction", .x), 
#                                               pattern = "summary", 
#                                               full.names = TRUE))

summary_interaction_files <- map(c("Visium", "sn_broad", "sn_fine_leaveOut_Br3974"), ~list.files(here("processed-data", "12_voomLmFit", "03_Clusterwise_voomLmFit_interaction", paste0("vlmf_",.x)), 
                                              pattern = "summary", 
                                              full.names = TRUE))

names(summary_interaction_files) <- datatypes

summary_interaction_tbl <- map2_dfr(summary_interaction_files, names(summary_interaction_files), 
                                function(files, name){map(files, ~read.csv(.x) |> 
                                                              mutate(datatype = name, 
                                                                     file = basename(.x))
                                                          )}
                                ) 

summary_tbl <- summary_tbl |> 
    filter(is.na(mod) | mod == "carrier") |> 
    bind_rows(summary_interaction_tbl)
    

sum_FDR05_tbl <- summary_tbl |> 
    group_by(Mod, datatype) |>
    summarize(sum_FDR05 = sum(n_FDR05))


sum_FDR05_tbl |>
    ggplot(aes(x = datatype, y = Mod)) +
    geom_tile(aes(fill = sum_FDR05)) +
    geom_text(aes(label = sum_FDR05))

