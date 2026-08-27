# snRNA-seq Processing, Clustering, and Annotation

Scripts are numbered in pipeline order. Each `.R` script generally has a matching `.sh` SLURM wrapper. Scripts that **save an `sce` object to disk** are marked with 💾 and noted below — these are the checkpoints later scripts build on.

## QC and preprocessing

* `00_sn_data_check` Prepare snRNA-seq data; build sample metadata from ancestry/phenotype info and collect Chromium run info.
* `01_get_droplet_scores` Run droplet detection (`DropletUtils::emptyDrops`) per sample to distinguish real nuclei from empty droplets.
* `02_droplet_QC` 💾 Read in droplet scores, filter out empty droplets, run doublet detection (`scDblFinder`) and QC metrics. Saves `processed-data/sce_objects/sce_postQC`.
* `03_GLM_PCA` 💾 Deviance-based feature selection and GLM-PCA dimensionality reduction on the post-QC SCE. Saves `processed-data/sce_objects/sce_uncorrected`.
* `04_plot_uncorrected_reduced_dims` Plot UMAP/t-SNE on the uncorrected (GLM-PCA) embedding, colored by QC metrics and marker genes.
* `05_harmony_correction` 💾 Batch-correct the GLM-PCA embedding with Harmony, then recompute UMAP/t-SNE. Saves `processed-data/spe_objects/sce_harmony.Rdata` and `processed-data/sce_objects/sce_harmony`.
* `06_plot_corrected_reduced_dims` Plot Harmony-corrected UMAP/t-SNE with categorical, continuous, and marker gene overlays.

## Preliminary clustering and cluster QC

* `07_cluster_sn_prelim` Build an SNN graph and run walktrap clustering on Harmony-corrected data to get preliminary clusters (long runtime, ~13h).
* `08_sctype_prelim` Annotate preliminary clusters with ScType (brain cell type database).
* `09_cluster_QC` Evaluate preliminary clusters using QC metrics (doublet score, mito %, library size, detected genes) and flag low-quality clusters/nuclei for removal.
* `10_reprocess_quality_sn` 💾 Drop low-quality clusters/nuclei, then re-normalize, re-run GLM-PCA and Harmony correction on the retained nuclei. Saves `processed-data/sce_objects/sce_reprocess`.
* `10.5_cluster_dotplot` Dotplots of literature marker genes across the preliminary broad cell type clusters.

## Final clustering and annotation

* `11_cluster_sn` Cluster the reprocessed, high-quality nuclei with SNN + walktrap across several values of *k* (e.g. 10, 15, 20).
* `12_cluster_eval` Evaluate clustering resolution across *k* via silhouette score to pick an optimal *k*.
* `13_sctype_final` Annotate the optimal (final) clustering with ScType, cell type database selectable via command-line argument.
* `15_cluster_update_sce` 💾 Integrate the optimal clustering + ScType annotations (and colors) into the SCE colData/metadata. Saves `processed-data/sce_objects/sce_ERC` — the main annotated snRNA-seq object used throughout the rest of the pipeline.

## Marker genes, pseudobulk, and spatial registration (broad cell types)

* `16_sn_MeanRatio` Compute mean-ratio marker gene statistics per cell type (`DeconvoBuddies::get_mean_ratio`).
* `17_sn_model_pseudobulk` Run `spatialLIBD::registration_wrapper` to fit cell type modeling results and generate pseudobulked data by cluster/sample.
* `18_sn_spatial_registration` Correlate snRNA-seq cell type modeling results against DLPFC and ERC spatial domains ("spatial registration").
* `19_sn_heatmaps` Heatmaps of pseudobulked expression across cell types/samples, including spatial/DLPFC enrichment correlations.
* `19.5_sn_dotplot` Dotplots of literature marker genes across broad cell types on the full `sce_ERC` object.
* `20_sn_pseudobulk` Aggregate cells across a chosen cluster variable (`scuttle::aggregateAcrossCells`) and CPM-normalize. Saves a pseudobulked **SingleCellExperiment** via `saveRDS()` (not HDF5) to `processed-data/04_snRNA-seq/20_sn_pseudobulk/{sce}_pseudobulk_only-{cluster_var}.rds` — this is pseudobulk output, not the full-resolution `sce`, but is flagged here since it is an `sce` object.
* `21_sn_pseudobulk_pca` PCA on the pseudobulked data; plot variance explained and PCs vs. covariates.
* `22_crumblr_sn_broad` Differential proportion analysis (Crumblr) on **broad** cell type proportions across samples/covariates.
* `22_crumblr_sn` Differential proportion analysis (Crumblr) on **fine-grained** cell type annotations, with phenotype covariates.
* `23_sn_model_exploration` Ad hoc exploration of modeling results, e.g. cell-type-specific DE for microglia subtypes.
* `24_external_data_check` Compare/register the snRNA-seq data against external reference datasets (e.g. Sestan EC) via `registration_wrapper`.

## Sub-clustering glia/subtypes

* `25_sn_cluster_subtype` Sub-cluster individual broad cell types (Astro, Endo, Micro, Oligo, OPC) with GLM-PCA + SNN clustering.
* `26_sn_subtype_check` Annotate sub-cluster results per cell type using marker expression, ScType, and manual inspection.
* `27_sn_glia_check` Compare glia re-clustering vs. sub-clustering results and consolidate glia subtype annotations (produces the annotation spreadsheet used downstream).
* `28_subcluster_update_sce` 💾 Add sub-cluster/glia annotations to the main SCE, update ancestry data, drop QC-failing clusters, and finalize `cell_type_anno`/`cell_type_broad`. Saves `processed-data/sce_objects/sce_ERC_subcluster` — the finalized, sub-clustered annotation object used by later scripts.
* `29_sn_subcluster_model_pseudobulk` Run `registration_wrapper` for modeling + pseudobulk generation on the sub-clustered cell types.
* `30_sn_subcluster_spatial_registration` Correlate sub-cluster modeling results with DLPFC layers, spatial ERC domains, and external references.
* `31_sn_subcluster_heatmap` Heatmaps of pseudobulked sub-cluster data with sample/cell-type annotations.
* `32_sn_subcluster_hierarchical_cluster` Hierarchical clustering of pseudobulked sub-cluster data to inspect clustering resolution/relationships.
* `33_sn_subcluster_summary` Summary tables and exploratory plots of the final sub-cluster annotations.
* `34_sn_subcluster_MeanRatio` Mean-ratio marker gene statistics for the sub-clustered cell types.
* `35_sn_subcluster_marker_modeling` Marker gene modeling and dotplots (e.g. AD risk genes) by cell type sub-type.
* `35.5_sn_subcluster_marker_modeling_OligoOPC` Same marker modeling, focused on Oligodendrocyte/OPC differentiation trajectories.
* `36_sn_subcluster_dotplot` Dotplots of sex-linked and risk genes across sub-cluster cell types and samples.
* `37_sn_subcluster_reducedDims` Compute/plot UMAP/t-SNE for individual sub-clustered cell types.
* `38_sn_subcluster_reducedDims_OligoOPC` Reduced dimensions + trajectory analysis (TSCAN) for combined Oligo + OPC populations.

## Other

* `00_lit_marker_genes` Compile and reconcile literature marker gene lists (e.g. PsychENCODE/AHBA PFC) for brain cell types, used by dotplot scripts throughout.

## SCE checkpoints at a glance

| Script | Saves | Path |
|---|---|---|
| `02_droplet_QC` | sce (HDF5) | `processed-data/sce_objects/sce_postQC` |
| `03_GLM_PCA` | sce (HDF5) | `processed-data/sce_objects/sce_uncorrected` |
| `05_harmony_correction` | sce (Rdata + HDF5) | `processed-data/spe_objects/sce_harmony.Rdata`, `processed-data/sce_objects/sce_harmony` |
| `10_reprocess_quality_sn` | sce (HDF5) | `processed-data/sce_objects/sce_reprocess` |
| `15_cluster_update_sce` | sce (HDF5) | `processed-data/sce_objects/sce_ERC` — main annotated object |
| `20_sn_pseudobulk` | sce (RDS, pseudobulked) | `processed-data/04_snRNA-seq/20_sn_pseudobulk/{sce}_pseudobulk_only-{cluster_var}.rds` |
| `28_subcluster_update_sce` | sce (HDF5) | `processed-data/sce_objects/sce_ERC_subcluster` — final sub-clustered object |
