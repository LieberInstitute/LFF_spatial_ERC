library(tidyverse)
library(here)
library(sessioninfo)

deg_path = here(
    'processed-data', '13_compile_DGE', '01_compile_DGE', 'Xenium',
    'DGE_results_carrier_Xenium_wSN.csv'
)
pair_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'unique_ligand_receptor_pairs.csv'
)
out_path = here(
    'processed-data', '24_xenium_liana', '07_DEG_overlap',
    'lr_deg_overlap.csv'
)

pair_df = read_csv(pair_path, show_col_types = FALSE)

#   Load DEG results and check whether any LR genes are significant in Xenium
deg = read_csv(deg_path, show_col_types = FALSE)
lr_genes = union(pair_df$ligand_complex, pair_df$receptor_complex)

#   No Xenium hits for any LR gene
stopifnot(
    deg |>
        filter(gene_name %in% lr_genes, vlmf_xenium_adj.P.Val < 0.05) |>
        nrow() == 0
)

#   For each LR gene, keep the cluster with the lowest SN adj p-value
deg_sn = deg |>
    filter(gene_name %in% lr_genes) |>
    select(cluster, gene_name, vlmf_sn_logFC, vlmf_sn_adj.P.Val) |>
    slice_min(vlmf_sn_adj.P.Val, by = gene_name, n = 1, with_ties = FALSE)

#   Join SN DEG info for ligand and receptor into a single row per LR pair
lr_deg = pair_df |>
    left_join(
        deg_sn |>
            select(
                gene_name,
                ligand_cluster = cluster,
                ligand_sn_logFC = vlmf_sn_logFC,
                ligand_sn_padj = vlmf_sn_adj.P.Val
            ),
        by = c('ligand_complex' = 'gene_name')
    ) |>
    left_join(
        deg_sn |>
            select(
                gene_name,
                receptor_cluster = cluster,
                receptor_sn_logFC = vlmf_sn_logFC,
                receptor_sn_padj = vlmf_sn_adj.P.Val
            ),
        by = c('receptor_complex' = 'gene_name')
    )

write_csv(lr_deg, out_path)

session_info()
