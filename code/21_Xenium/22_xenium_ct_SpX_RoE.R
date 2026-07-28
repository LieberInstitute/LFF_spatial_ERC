## Louise Huuki-Myers, July 2026
## Check cell types in Xenium data for divergence from Ro/e

#### Set Up ####

library("qs2")
library("SpatialExperiment")
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("ComplexHeatmap")
library("circlize")



plot_dir <- here("plots", "21_Xenium", "22_xenium_ct_SpX_RoE")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "21_Xenium", "22_xenium_ct_SpX_RoE")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load data ####

message(Sys.time(), " - Load SPE data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))
# filter to singlet cells
spe <- spe[,spe$spot_class == "singlet"]

table(spe$cell_type_anno == "Oligo.3")

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
SpX_colors <- metadata(spe)$SpX_colors

cell_v_SpX <- table(SpX = spe$SpX, cell_type_anno = spe$cell_type_anno)

#### Ro/e enrichment ####
## Ratio of observed to expected cell frequency (Zhang et al. 2018 Nature;
## Guo et al. 2018 Nat Med). Ro/e = 1 means a cell type is distributed across
## SpX domains exactly in proportion to each domain's overall size (i.e. no
## spatial preference); Ro/e > 1 = enriched, Ro/e < 1 = depleted.
##
## chisq.test()$expected conveniently returns exactly this null expectation
## (row_total * col_total / grand_total), so we reuse it rather than
## calculating domain size shares by hand.

message(Sys.time(), " - Compute Ro/e for cell_type_anno x SpX")
chisq_result <- chisq.test(cell_v_SpX)
# NOTE: chisq.test may warn "approximation may be incorrect" for rare cell
# types with very low expected counts (e.g. OPC.1/OPC.3) - this only affects
# the overall test statistic, not the expected counts used below.

RoE_cutoff <- 1.5 # fold-change threshold for calling enrichment/depletion

cell_v_SpX_RoE <- as.data.frame(cell_v_SpX) |>
  as_tibble() |>
  rename(n_cell = Freq) |>
  mutate(
    expected_n = as.vector(chisq_result$expected),
    std_resid = as.vector(chisq_result$stdres),
    # Haldane-Anscombe-style pseudocount avoids 0/0 and log2(0) = -Inf when a
    # cell type has zero cells in a given domain (matches the per-donor calc below)
    RoE = (n_cell + 0.5) / (expected_n + 0.5),
    log2_RoE = log2(RoE),
    enrichment = case_when(
      RoE >= RoE_cutoff ~ "enriched",
      RoE <= 1 / RoE_cutoff ~ "depleted",
      TRUE ~ "neutral"
    )
  ) |>
  group_by(cell_type_anno) |>
  mutate(prop_SpX = n_cell / sum(n_cell)) |> # density: this cell type's cells, split across domains (sums to 1 per cell type)
  group_by(SpX) |>
  mutate(prop_cell_type = n_cell / sum(n_cell)) |> # composition: this domain's cells, split across cell types (sums to 1 per domain)
  ungroup() |>
  arrange(cell_type_anno, desc(RoE))

message(
  Sys.time(), " - N cell type x SpX pairs enriched (RoE >= ", RoE_cutoff, "): ",
  sum(cell_v_SpX_RoE$enrichment == "enriched")
)

## IMPORTANT CAVEAT: std_resid/chisq p-values assume independent cells, which
## is violated here since cells are nested within donors (pseudoreplication).
## Treat RoE as a descriptive effect size and std_resid as a rough guide only
## - the by-donor test below is the calibrated significance measure.

#### Per-donor Ro/e (accounts for donor pseudoreplication) ####
## Same Ro/e logic, computed separately within each donor's own 2-way table
## (spe$sample_id), so that donor - not cell - is the unit of replication for
## a formal significance test.

message(Sys.time(), " - Compute per-donor Ro/e")

donor_tab <- table(
  sample_id = spe$sample_id,
  SpX = spe$SpX,
  cell_type_anno = spe$cell_type_anno
)

min_n_cell_donor <- 100 # skip donors with too few captured cells to trust their expected counts

donor_RoE <- map_dfr(dimnames(donor_tab)$sample_id, function(id) {
  tab_d <- donor_tab[id, , ]
  if (sum(tab_d) < min_n_cell_donor) return(NULL)
  expected_d <- chisq.test(tab_d)$expected
  as.data.frame(tab_d) |>
    as_tibble() |>
    rename(n_cell = Freq) |>
    mutate(
      sample_id = id,
      expected_n = as.vector(expected_d),
      # Haldane-Anscombe-style pseudocount avoids 0/0 and Inf when a cell
      # type is entirely absent from a domain (or donor) by chance
      RoE = (n_cell + 0.5) / (expected_n + 0.5),
      log2_RoE = log2(RoE)
    )
})

#### One-sample test per cell type x SpX, across donors ####
## H0: mean log2(RoE) across donors = 0 (no consistent enrichment/depletion).
## Pairs with fewer than min_n_donor contributing donors are left untested
## (NA) rather than given an unreliable p-value from too few observations.

min_n_donor <- 5

donor_RoE_test <- donor_RoE |>
  group_by(cell_type_anno, SpX) |>
  summarize(
    n_donor = n(),
    mean_log2_RoE = mean(log2_RoE),
    p_value = if (n_donor >= min_n_donor) wilcox.test(log2_RoE, mu = 0)$p.value else NA_real_,
    .groups = "drop"
  ) |>
  mutate(padj = p.adjust(p_value, method = "BH"))

message(
  Sys.time(), " - N cell type x SpX pairs with padj < 0.05 (by-donor test): ",
  sum(donor_RoE_test$padj < 0.05, na.rm = TRUE)
)

## join by-donor test results back into the pooled table
cell_v_SpX_RoE <- cell_v_SpX_RoE |>
  left_join(donor_RoE_test, by = c("cell_type_anno", "SpX")) |>
  mutate(sig_by_donor = !is.na(padj) & padj < 0.05)

#### Primary domain per cell type ####
## One row per cell_type_anno: the SpX domain with the strongest Ro/e

primary_domain_tab <- cell_v_SpX_RoE |>
  group_by(cell_type_anno) |>
  slice_max(RoE, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(cell_type_anno, primary_SpX = SpX, RoE, prop_SpX, enrichment, n_donor, padj, sig_by_donor)

#### Save results ####

write_csv(cell_v_SpX_RoE, file.path(data_dir, "cell_v_SpX_RoE.csv"))
write_csv(primary_domain_tab, file.path(data_dir, "cell_type_primary_SpX.csv"))
qs_save(cell_v_SpX_RoE, file.path(data_dir, "cell_v_SpX_RoE.qs2"))

#### Heatmap of cell type x SpX: counts, density, and Ro/e enrichment ####
## Row/column annotations (SpX colors, fine cell_type colors) and the broad
## cell type column split (stripping ".N" suffix) are shared across all
## three heatmaps so they line up directly for comparison.

stopifnot(all(colnames(cell_v_SpX) %in% names(cell_type_colors$anno)))
stopifnot(all(rownames(cell_v_SpX) %in% names(SpX_colors)))

## density: proportion of each cell type's cells across SpX domains (columns sum to 1)
cell_v_SpX_prop <- prop.table(cell_v_SpX, margin = 2)

## long -> wide Ro/e matrix, ordered to match the original cell_v_SpX dimnames
## so annotation colors line up correctly
log2_RoE_mat <- cell_v_SpX_RoE |>
  select(SpX, cell_type_anno, log2_RoE) |>
  pivot_wider(names_from = cell_type_anno, values_from = log2_RoE) |>
  column_to_rownames("SpX") |>
  as.matrix()
log2_RoE_mat <- log2_RoE_mat[rownames(cell_v_SpX), colnames(cell_v_SpX)]

## asterisk = enriched (pooled RoE >= cutoff) AND significant across donors
sig_mat <- cell_v_SpX_RoE |>
  mutate(sig_mark = if_else(enrichment == "enriched" & sig_by_donor, "*", "")) |>
  select(SpX, cell_type_anno, sig_mark) |>
  pivot_wider(names_from = cell_type_anno, values_from = sig_mark) |>
  column_to_rownames("SpX") |>
  as.matrix()
sig_mat <- sig_mat[rownames(cell_v_SpX), colnames(cell_v_SpX)]

## diverging color scale for Ro/e, centered at log2(RoE) = 0
RoE_max <- max(abs(log2_RoE_mat), na.rm = TRUE)
RoE_col_fun <- circlize::colorRamp2(c(-RoE_max, 0, RoE_max), c("steelblue", "white", "firebrick"))

## create annotations (shared across all three heatmaps - dimensions match)
SpX_row_ha <- rowAnnotation(
  SpX = rownames(cell_v_SpX),
  col = list(SpX = SpX_colors),
  show_legend = FALSE
)

cell_type_col_ha <- HeatmapAnnotation(
  cell_type = colnames(cell_v_SpX),
  col = list(cell_type = cell_type_colors$anno),
  show_legend = FALSE
)

sig_legend <- ComplexHeatmap::Legend(
  labels = "enriched (RoE >= cutoff) &\nsignificant by-donor (padj < 0.05)",
  type = "points", pch = "*", legend_gp = grid::gpar(col = "black")
)

## PLOT HEATMAPS
pdf(here(plot_dir, "Xenium_bansky_SpX_v_cell_type_heatmap_broad.pdf"), width = 10)

## raw counts
ComplexHeatmap::Heatmap(cell_v_SpX,
  name = "n singlet cells",
  col = c("black", viridisLite::plasma(100)),
  column_split = gsub("\\..*", "", colnames(cell_v_SpX)),
  cluster_rows = FALSE,
  left_annotation = SpX_row_ha,
  bottom_annotation = cell_type_col_ha
)

## density (proportion SpX, cell type cols sum to 1)
ComplexHeatmap::Heatmap(cell_v_SpX_prop,
  name = "prop SpX\nsinglet cells",
  col = c("black", viridisLite::plasma(100)),
  column_split = gsub("\\..*", "", colnames(cell_v_SpX_prop)),
  cluster_rows = FALSE,
  left_annotation = SpX_row_ha,
  bottom_annotation = cell_type_col_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sig_mat[i, j], x, y, gp = grid::gpar(fontsize = 10))
  }
)

## Ro/e enrichment, with by-donor significance stars
RoE_ht <- ComplexHeatmap::Heatmap(log2_RoE_mat,
  name = "log2(Ro/e)",
  col = RoE_col_fun,
  column_split = gsub("\\..*", "", colnames(log2_RoE_mat)),
  cluster_rows = FALSE,
  left_annotation = SpX_row_ha,
  bottom_annotation = cell_type_col_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(sig_mat[i, j], x, y, gp = grid::gpar(fontsize = 10))
  }
)
# ComplexHeatmap::draw(RoE_ht, annotation_legend_list = list(sig_legend))
ComplexHeatmap::draw(RoE_ht)

dev.off()


