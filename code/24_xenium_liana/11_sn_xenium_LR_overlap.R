library(tidyverse)
library(here)
library(sessioninfo)

sn_path = here(
    'processed-data', '17_cell_cell_comm', '10_LR_DEG_check',
    'Liana_LR_DEGs_sn_fine.csv'
)
xenium_path = here(
    'processed-data', '24_xenium_liana', '05_top_scatter',
    'pair_level_comparison.csv'
)
out_dir = here('processed-data', '24_xenium_liana', '11_sn_xenium_LR_overlap')
out_path = file.path(out_dir, 'sn_xenium_LR_overlap.csv')

p_cut = 0.1
top_n = 5
key_cols = c('source', 'target', 'ligand_complex', 'receptor_complex')

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#   The sn table only contains LR pairs with a DEG on either side, so a pair
#   missing from it means "not DEG-linked" among the sn-tested pairs
sn_df = read_csv(sn_path, show_col_types = FALSE) |>
    select(all_of(key_cols), L_DEG, R_DEG) |>
    mutate(sn_tested = TRUE)

xenium_df = read_csv(xenium_path, show_col_types = FALSE) ## only on jhpce

stopifnot(
    all(c(key_cols, 'is_nom_sig', 'lr_mean_diff', 'abs_lr_mean_diff') %in%
        colnames(xenium_df))
)

xenium_df = xenium_df |>
    mutate(xenium_tested = TRUE)

overlap_df = full_join(xenium_df, sn_df, by = key_cols) |>
    group_by(ligand_complex, receptor_complex, 
             source_O3 = source == "Oligo.3", 
             target_O3 = target == "Oligo.3") |>
    mutate(
        sn_tested = replace_na(sn_tested, FALSE),
        xenium_tested = replace_na(xenium_tested, FALSE),
        ligand_sn_DEG = replace_na(L_DEG, FALSE),
        receptor_sn_DEG = replace_na(R_DEG, FALSE),
        sn_DEG_support = ligand_sn_DEG | receptor_sn_DEG,
        #   Rank by absolute percent difference among Xenium-tested pairs
        xen_diff_rank = if_else(
            xenium_tested,
            min_rank(desc(abs_lr_mean_diff)),
            NA_integer_
        ),
        xen_ttest_support = replace_na(is_nom_sig, FALSE),
        xen_top_diff_support = replace_na(xen_diff_rank <= top_n, FALSE),
        n_tests_supported = sn_DEG_support + xen_ttest_support +
            xen_top_diff_support
    ) |>
    select(
        sender = source, ligand = ligand_complex, target,
        receptor = receptor_complex, n_tests_supported,
        sn_DEG_support, ligand_sn_DEG, receptor_sn_DEG, sn_tested,
        xen_ttest_support, is_nom_sig, xenium_tested,
        xen_top_diff_support, xen_diff_rank, lr_mean_diff,
        abs_lr_mean_diff, `lr_mean_E2+`, `lr_mean_E4+`, higher_in
    ) |>
    ungroup() |>
    arrange(desc(n_tests_supported), xen_diff_rank) 

overlap_df |> ungroup() |> count(n_tests_supported)
overlap_df |> ungroup() |> count(xen_top_diff_support)
overlap_df |> filter(receptor == "ERBB3", ligand == "NRG3", xen_top_diff_support)

overlap_df |> filter(n_tests_supported > 1)

overlap_df |> 
    ungroup() |>
    filter(receptor == "ERBB3", target == "Oligo.3", ligand == "NRG3") |>
    # filter(receptor == "NRXN1", target == "Oligo.3", ligand == "NLGN1") |>
    select(sender, ligand, target, receptor, sn_DEG_support, xen_ttest_support, xen_top_diff_support, xen_diff_rank) |>
    print(n = 34)

overlap_df |> 
    filter(sn_DEG_support) |>
    select(sender, ligand, target, receptor, sn_DEG_support, xen_ttest_support, xen_top_diff_support, xen_diff_rank) |>
    print(n = 20)

write_csv(overlap_df, out_path)

overlap_df |> count(n_tests_supported)

session_info()
