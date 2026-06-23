#' Find the nearest neighbor of a certain cell type for the reference cell
#'
#' @param spe 
#' @param reference_barcode 
#' @param neighbor 
#'
#' @returns
#' @export
#'
#' @examples
#' message(Sys.time(), " - Load SPE data")
#' spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))
#' # filter to singlet cells
#' spe <- spe[,spe$spot_class == "singlet"]
#' 
#' ## filter to one sample to test 
#' spe_test <- spe[,spe$BrNum == "Br1039"]
#' 
#' ncol(spe_test) #23k
#' 
#' nearest_cell_type_neighbor(spe = spe_test,
#'                            reference_barcode = "Br1039_aabdmomo-1",
#'                            neighbor = "Astro",
#'                            neighbor_ct_col = "cell_type_broad")
#'
nearest_cell_type_neighbor <- function(spe, 
                                       reference_barcode, 
                                       neighbor = "Astro", 
                                       neighbor_ct_col = "cell_type_broad") {
    
    # extract coords and cell type
    coords <- SpatialExperiment::spatialCoords(spe)
    cell_type <- colData(spe)[[neighbor_ct_col]]   # adjust column name if needed
    barcodes  <- spe$cell_id ## use cell_id for donor specific barcode
    
    # get query coords for the reference barcode
    ref_idx <- match(reference_barcode, barcodes)
    if (is.na(ref_idx)) stop("reference_barcode not found in spe")
    ref_coords <- coords[ref_idx, , drop=FALSE]
    
    # get candidate coords for the neighbor cell type
    neighbor_idx <- which(cell_type == neighbor)
    if (length(neighbor_idx) == 0) stop("No cells found for neighbor cell type: ", neighbor)
    neighbor_coords <- coords[neighbor_idx, , drop=FALSE]
    
    # nearest neighbor search (k=1)
    nn <- RANN::nn2(data=neighbor_coords, query=ref_coords, k=1)
    
    # resolve back to barcode
    nearest_barcode <- barcodes[neighbor_idx[nn$nn.idx[1]]]
    nearest_dist <- nn$nn.dists[1]
    
    list(reference_barcode = reference_barcode, neightbor_barcode=nearest_barcode, distance=nearest_dist)
}


#' Find the nearest neighbors for all refrences from a sample
#'
#' @param spe 
#' @param donor 
#' @param reference 
#' @param reference_ct_col 
#' @param neighbor 
#' @param neighbor_ct_col 
#'
#' @returns
#' @export
#'
#' @examples
#' test <- map_nearest_cell_type_neighbor(spe,
#'                                        donor = "Br1039",
#'                                        reference = "Oligo.3",
#'                                        reference_ct_col = "cell_type_anno",
#'                                        neighbor = "Astro",
#'                                        neighbor_ct_col = "cell_type_broad")
#'                                        
map_nearest_cell_type_neighbor <- function(spe, 
                                           donor = "Br1039",
                                           reference = "Oligo.3",
                                           reference_ct_col = "cell_type_anno",
                                           neighbor = "Astro", 
                                           neighbor_ct_col = "cell_type_broad"){
    
    ## filter to one donor
    spe <- spe[,spe$BrNum == donor]
    
    spe_ref <- spe[,spe[[reference_ct_col]] == reference]
    # message(Sys.time(), sprintf(" - %s reference cells: %s, ncells %i",donor, reference, ncol(spe_ref)))
    
    neighbor_df <- map_dfr(spe_ref$cell_id, ~nearest_cell_type_neighbor(spe = spe, 
                                                         reference_barcode = .x, 
                                                         neighbor = neighbor, 
                                                         neighbor_ct_col = neighbor_ct_col)) |>
        mutate(BrNum = donor, .before=1)
    
    
    return(neighbor_df)
    
}

#' Find the nearest neighbors for all references across all donors in parallel
#'
#' @param spe A SpatialExperiment object containing all donors
#' @param donors Character vector of donor IDs (BrNum values). If NULL, uses all unique BrNum values.
#' @param reference Cell type label for the reference cells
#' @param reference_ct_col Column name in colData for reference cell type
#' @param neighbor Cell type label for the neighbor cells
#' @param neighbor_ct_col Column name in colData for neighbor cell type
#' @param BPPARAM A BiocParallelParam object. Defaults to bpparam().
#'
#' @returns A data.frame with columns: BrNum, reference_barcode, neighbor_barcode, distance
#' @export
#'
#' @examples
#' library(BiocParallel)
#'
#' ## serial (for testing)
# neighbor_df <- bpmap_nearest_cell_type_neighbor(spe,
#                                                 donors = c("Br1039", "Br1556"),
#                                                 reference = "Oligo.3",
#                                                 reference_ct_col = "cell_type_anno",
#                                                 neighbor = "Astro",
#                                                 neighbor_ct_col = "cell_type_broad",
#                                                 BPPARAM = BiocParallel::SerialParam())
#'
#' ## parallel across donors
#' neighbor_df <- bpmap_nearest_cell_type_neighbor(spe,
#'                                                 donors = c("Br1039", "Br1556", "Br1706", "Br2582"),
#'                                                 reference = "Oligo.3",
#'                                                 reference_ct_col = "cell_type_anno",
#'                                                 neighbor = "Astro",
#'                                                 neighbor_ct_col = "cell_type_broad",
#'                                                 BPPARAM = BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE))
#'
bpmap_nearest_cell_type_neighbor <- function(spe,
                                             donors = NULL,
                                             reference = "Oligo.3",
                                             reference_ct_col = "cell_type_anno",
                                             neighbor = "Astro",
                                             neighbor_ct_col = "cell_type_broad",
                                             BPPARAM = bpparam()) {
    
    if (is.null(donors)) donors <- unique(spe$BrNum)
    message(Sys.time(), sprintf(" - Running across %i donors: %s", length(donors)))
    
    results <- BiocParallel::bplapply(
        donors,
        FUN = map_nearest_cell_type_neighbor,
        spe = spe,
        reference = reference,
        reference_ct_col = reference_ct_col,
        neighbor = neighbor,
        neighbor_ct_col = neighbor_ct_col,
        BPPARAM = BPPARAM
    )
    
    dplyr::bind_rows(results)
}
