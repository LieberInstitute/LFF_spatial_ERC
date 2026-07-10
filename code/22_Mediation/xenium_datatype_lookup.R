## pb_fn / DE_fn resolution + caching, shared across the mediation
## scripts (03_Mediation_Xenium_pairs.R, 04_Mediation_Xenium_plots.R).
## Pulled out into its own file so the two scripts can't silently drift
## apart on which file a given datatype key points to.
##
## Provides:
##   pb_lookup / get_pb_fn() / load_pb_cached()  -- pseudobulk objects
##   DE_lookup / get_DE_fn() / load_DE_cached()  -- compiled DE results
##     (used as mediator_stats for plot_mediator_panel=TRUE)
##
## Both lookups use the SAME datatype_key values, so a given
## mediator_datatype/outcome_datatype resolves consistently to both its
## pseudobulk object and its DE_data table.

pb_dir <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep")
pb_lookup <- tribble(
    ~datatype_key,                               ~pb_fn,
    "Xenium_cell_type_anno",                     file.path(pb_dir, "spe_xenium_pseudo_DGE-cell_type_anno.RDS"),
    "Xenium_cell_type_anno_SpX",                 file.path(pb_dir, "spe_xenium_pseudo_DGE-cell_type_anno_SpX.RDS"),
    "Xenium_Oligo.3_Astro",                      file.path(pb_dir, "spe_xenium_pseudo_DGE-Oligo.3_Astro.RDS"),
    "Xenium_Oligo.3_Astro_SpX",                  file.path(pb_dir, "spe_xenium_pseudo_DGE-Oligo.3_Astro_SpX.RDS")
)

get_pb_fn <- function(datatype_key) {
    hit <- pb_lookup$pb_fn[pb_lookup$datatype_key == datatype_key]
    if (length(hit) != 1) stop(sprintf("Unrecognized datatype key '%s' -- add it to pb_lookup.", datatype_key))
    hit
}

.pb_cache <- new.env()
load_pb_cached <- function(fn) {
    if (!exists(fn, envir = .pb_cache, inherits = FALSE)) {
        message(Sys.time(), " - loading pseudobulk: ", basename(fn))
        sce <- readRDS(fn)
        sce$APOE_carrier_syn <- gsub("\\+", "", sce$APOE_carrier)
        assign(fn, sce, envir = .pb_cache)
    }
    get(fn, envir = .pb_cache, inherits = FALSE)
}

## mediator_stats (compiled DGE results) resolution + caching, mirroring
## pb_lookup/load_pb_cached -- these are the DE_data tables used for the
## mediator panel (plot_mediator_panel=TRUE), and there's one per
## datatype, same as the pseudobulk objects. Using the wrong one (e.g.
## always the plain cell_type_anno table) silently fails to find rows for
## the APOE/spatially-stratified mediator scenarios.
DE_dir <- here("processed-data", "13_compile_DGE", "01_compile_DGE")
DE_lookup <- tribble(
    ~datatype_key,                 ~DE_fn,
    "Xenium_cell_type_anno",       file.path(DE_dir, "Xenium_cell_type_anno", "DGE_results_carrier_Xenium_cell_type_anno.Rds"),
    "Xenium_cell_type_anno_SpX",   file.path(DE_dir, "Xenium_SpX", "DGE_results_carrier_Xenium_SpX.Rds"),
    "Xenium_Oligo.3_Astro",        file.path(DE_dir, "Xenium_Oligo.3_Astro", "DGE_results_carrier_Xenium_Oligo.3_Astro.Rds"),
    "Xenium_Oligo.3_Astro_SpX",    file.path(DE_dir, "Xenium_Oligo.3_Astro_SpX", "DGE_results_carrier_Xenium_Oligo.3_Astro_SpX.Rds")
)

get_DE_fn <- function(datatype_key) {
    hit <- DE_lookup$DE_fn[DE_lookup$datatype_key == datatype_key]
    if (length(hit) != 1) stop(sprintf("Unrecognized datatype key '%s' -- add it to DE_lookup.", datatype_key))
    hit
}

.DE_cache <- new.env()
load_DE_cached <- function(fn) {
    if (!exists(fn, envir = .DE_cache, inherits = FALSE)) {
        message(Sys.time(), " - loading DE_data: ", basename(fn))
        assign(fn, readRDS(fn), envir = .DE_cache)
    }
    get(fn, envir = .DE_cache, inherits = FALSE)
}
