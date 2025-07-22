# LFF_spatial_ERC
LFF_spatial_ERC

## Internal

JHPCE location: /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/


## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔭 [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites to explore the ERC data:

- 🔭 [LFF_ERC_Visium](https://libd.shinyapps.io/LFF_ERC_Visium/):
    Visium data for 31 samples post-QC w/ clustering and modeling results
    
- 🔭 [LFF_ERC_Visium_QC](https://libd.shinyapps.io/LFF_ERC_Visium_QC/):
    Visium data for 31 samples pre-QC drops (including out-tissue spots)
    
- 👀 [iSEE](https://libd.shinyapps.io/LFF_ERC_snRNA-seq/):
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

**FASTQ**

`raw-data/FASTQ_snRNAseq`

### Visium

**SpatialExperiment**

project path: `"processed-data/spe_objects/spe_ERC_annotated"`
JHPCE path: `"/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/spe_objects/spe_ERC_annotated"`

Not on github due to size (4G)

```
## Load HD5F spe
spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
# class: SpatialExperiment 
# dim: 30494 122202 

```

**FASTQ** 

`raw-data/FASTQ`

