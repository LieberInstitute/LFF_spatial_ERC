# LFF_spatial_ERC
LFF_spatial_ERC

## Internal

JHPCE location: `/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/`


## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔭 [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites to explore the ERC data:

- 🔭 [LFF_ERC_Visium](https://interactive.libd.org/LFF_ERC_spatialLIBD-app/):
    Visium data for 31 samples post-QC w/ clustering and modeling results
    
- 🔭 [LFF_ERC_Visium_QC](https://libd.shinyapps.io/LFF_ERC_Visium_QC/):
    Visium data for 31 samples pre-QC drops (including out-tissue spots)
    
- 👀 [iSEE](https://interactive.libd.org/LFF_ERC_iSEE-app/):
    snRNA-seq data for 31 samples
    
## Main Data

### snRNA-seq

**SingleCellExperiment** 

project path: `"processed-data/sce_objects/sce_ERC_subcluster"`

JHPCE path: `"/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/sce_objects/sce_ERC_subcluster"`

Not on github due to size (36G)

```
## Load HD5F sce
sce <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "sce_objects", "sce_ERC_subcluster"))
# class: SingleCellExperiment 
# dim: 38606 122004

```

**Pseudobulked snRNA-seq**

cell type fine: `processed-data/08_pseudoBulkDGE_sn/01_pseudobulk_data_sn/sce_pseudo_DGE-cell_type_anno.RDS`
cell type broad: `processed-data/08_pseudoBulkDGE_sn/01_pseudobulk_data_sn/sce_pseudo_DGE-cell_type_broad.RDS`

**FASTQ**

`raw-data/FASTQ_snRNAseq`

**Sample ID table**
`processed-data/04_snRNA-seq/erc_sn_sample_info.csv`

### Visium

**SpatialExperiment**

project path: `"processed-data/spe_objects/spe_ERC_annotated"`

JHPCE path: `/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/spe_objects/spe_ERC_annotated`

pre-QC version: `processed-data/spe_objects/spe_raw.rds`

Not on github due to size (4G)

```
## Load HD5F spe
spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
# class: SpatialExperiment 
# dim: 30494 122202 

```

**Pseudobulked Visium data**

`processed-data/09_pseudoBulkDGE_Visium/01_pseudobulk_data_Visium/spe_pseudo_DGE.RDS`

**FASTQ** 

`raw-data/FASTQ`

**Images**

`raw-data/FASTQ/Images`

**SAMPLE ID table**

`processed-data/02_build_spe/sample_info.csv`

