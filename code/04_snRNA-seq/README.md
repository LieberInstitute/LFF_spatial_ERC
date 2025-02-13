
# Quality Control and Clustering for snRNA-seq data

* `00_sn_data_check` Prep snRNA-seq data, Build metadata and collect sample info

* `01_get_droplet_scores` Run droplet detection (DropletUtils::emptyDrops) by sample

* `02_droplet_QC` Read in output from 01_get_droplet_scores, filter out empty droplets [hdf5_sce_postQC]

* `03_GLM_PCA` Read in postQC SCE created in 02_droplet_qc, Run GLM PCA [sce_uncorrected]

* `04_plot_uncorrected_reduced_dims` Plot uncorrected reduced dimensions for snRNA-seq data

* `05_harmony_correction` Batch correction with Harmony [sce_harmony]

* `06_plot_corrected_reduced_dims`  Plot Harmony corrected reduced dimensions for snRNA-seq data

* `07_sn_prelim_cluster` Use single nearest neighbors + walk trap to cluster single nuc data k=10 (Long runtime ~13h)

* `08_sctype_prelim` Annotate preliminary clusters with SC Type

* `09_cluster_qc` <- `11_cluster_QC.R`

CLEAN ^^^

* `10_GLM_PCA_Harmony2` Read in cluster QC SCE created in Run GLM PCA

* `11_cluster_sn` Cluster high quality nuclei over many values of k

* `12_cluster_eval` Evaluate silhouette score for clustering

* `13_sctype_fine` Annotate optimal clustering 

* `14_finalize_snRNAseq` ?




13_sn_MeanRatio.R
14_sn_model_pseudobulk.R
15_cluster_update_sce.R
15_sn_spatial_registration.R