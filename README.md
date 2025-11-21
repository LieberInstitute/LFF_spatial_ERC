# LFF Spatial ERC

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17108407.svg)](https://doi.org/10.5281/zenodo.17108407)

## Overview

Welcome to the `LFF_spatial_ERC` project! 🧠

In this work we built a transcriptional profile of the human entorhinal
cortex (ERC) with single nucleus RNA sequencing (snRNA-seq) and spatially
resolved transcriptomics (SRT) from 31 neurotypical post-mortem donors with 
diverse risk for Alzheimer's disease (AD). SRT data was generated with 10x 
Genomics Visium, snRNA-seq was generated with 10x Genomics Chromium. 

We then assessed changes in gene expression associated with risk of AD by 
APOE carrier status (E2+ vs. E4+), as well as *APOE* changes specific to ancestry 
groups (African ancestry or European ancestry), and sex. 

This work is being was performed by members of [Leonardo
Collado-Torres](http://lcolladotor.github.io/), [Kristen
Maynard](https://www.libd.org/team/kristen-maynard-phd/), and [Keri
Martinowich](https://www.libd.org/team/keri-martinowich-phd/), teams at the
[Lieber Institute for Brain Development](libd.org) as well as Mina Ryten's group
from [UK DRI Cambridge](https://rytenlab.com/RytenLab/About). 

Thank you for your interest in our work!

## Study Design
<img src="http://research.libd.org/LFF_spatial_ERC/img/ERC experiment V2.png" width="1000px" align="left" />

Schematic depicting cohort composition stratified across *APOE* 
genotype, sex, and ancestry as well as experimental design for generation of 
snRNA-seq and SRT data from postmortem human ERC tissue. 

## Interactive Websites

We provide the following interactive websites to explore the ERC data:

- 🔭 [LFF_ERC_Visium](https://interactive.libd.org/LFF_ERC_spatialLIBD-app/):
    Visium data for 31 samples post-QC w/ clustering and modeling results
    
- 🔭 [LFF_ERC_Visium_QC](https://libd.shinyapps.io/LFF_ERC_Visium_QC/):
    Visium data for 31 samples pre-QC drops (including out-tissue spots)
    
- 👀 [iSEE](https://interactive.libd.org/LFF_ERC_iSEE-app/):
    snRNA-seq data for 31 samples

All of these interactive websites are powered by open source software,
namely:

- 🔭 [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) See 
[spatialLIBD website](https://research.libd.org/spatialLIBD/) for more information
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

## Citing our work

Please cite our pre-print [10.1101/2025.11.20.689483](https://doi.org/10.1101/2025.11.20.689483)
if you use data from this project.

> Huuki-Myers, L. A., Divecha, H. R., Bach, S. V., Valentine, M. R., Eagles, N. J.,
> Mulvey, B., Bharadwaj, R. A., Zhang, R., Evans, J. R., Grant-Peters, M., Miller,
> R. A., Kleinman, J. E., Han, S., Hyde, T. M., Page, S. C., Weinberger, D. R., 
> Martinowich, K., Ryten, M., Maynard, K. R. & Collado-Torres, L. 
> apoe E4 alzheimer’s risk converges on an oligodendrocyte subtype in the human 
> entorhinal cortex. BioRxiv (2025). doi:10.1101/2025.11.20.689483

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.


    @article{Huuki-Myers_Divecha_Bach_Valentine_Eagles_Mulvey_Bharadwaj_Zhang_Evans_Grant-Peters_et al._2025, 
    title={apoe E4 alzheimer’s risk converges on an oligodendrocyte subtype in the human entorhinal cortex}, 
    DOI={10.1101/2025.11.20.689483}, 
    abstractNote={<p>The entorhinal cortex (ERC) is implicated in early 
        progression of Alzheimer’s disease (AD). Here we investigated the impact of 
        established biological risk factors for AD, including APOE genotype 
        (E2 versus E4 alleles), sex, and ancestry, on gene expression in the human 
        ERC. We generated paired spatially-resolved transcriptomics (SRT) and 
        single-nucleus RNA sequencing data (snRNA-seq) in postmortem human ERC 
        tissue from middle aged brain donors with no history of AD. APOE-dependent 
        changes in gene expression predominantly mapped to a 
        transcriptionally-defined oligodendrocyte subtype, which varied 
        substantially with ancestry, and suggested differences in oligodendrocyte 
        differentiation and myelination. Integration of SRT and snRNA-seq data 
        identified a common gene expression signature associated with APOE genotype,
        which we localized to the same oligodendrocyte subtype and a white matter 
        spatial domain. This suggests that AD risk in ERC may be associated with 
        disrupted oligodendrocyte function, potentially contributing to future 
        neurodegeneration.</p>}, 
    journal={BioRxiv}, 
    author={Huuki-Myers, Louise A. and Divecha, Heena R. and Bach, Svitlana V. 
    and Valentine, Madeline R. and Eagles, Nicholas J. and Mulvey, Bernard and 
    Bharadwaj, Rahul A. and Zhang, Ruth and Evans, James R. and Grant-Peters, 
    Melissa and et al.}, 
    year={2025}, 
    month={Nov}}
    

## Data & Code avalibility 

Data and code for this project are available on [github](https://github.com/LieberInstitute/LFF_spatial_ERC).

Organization of code, data, and plots follows our team's [project template](https://github.com/LieberInstitute/template_project).

Raw data generated as a part of this study have been deposited on the Gene 
Expression Omnibus (GEO) with accession numbers GSE307990 [(SRT data)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307990) 
and GSE308007 [(snRNA-seq data)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GGSE308007).

Processed data files (R objects) used to make the interactive apps (logcounts only) 
can be found on our [Globus Endpoints](https://research.libd.org/globus/) under the 
heading `jhpce#LFF_ERC`. 

R objects with both the counts and logcounts can be 
downloaded through the fetch_data() function from [spatialLIBD website](https://research.libd.org/spatialLIBD/) 
v1.23.1 or newer. 

Source code is under the MIT License and permanently archived on 
[Zenodo](https://zenodo.org/records/17108408). The GitHub repository includes 
log files with specific software versions used for each analysis.

## Funding

This project was supported thanks to the Carol and Gene Ludwig Award for 
Neurodegeneration Research (DRW) from the [Ludwig Famiy Foundation](https://www.ludwigfamilyfoundation.org/) 
and the [Lieber Institute for Brain Development](https://www.libd.org/).


## Main Data

### snRNA-seq

**SingleCellExperiment** 

* project path: `"processed-data/sce_objects/sce_ERC_subcluster"` (Not on github due to large file size: 36G)

```
## Load HD5F sce
sce <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "sce_objects", "sce_ERC_subcluster"))
# class: SingleCellExperiment 
# dim: 38606 122004
```

** Input files **

* FASTQ: `raw-data/FASTQ_snRNAseq`

* Sample ID table: `processed-data/04_snRNA-seq/erc_sn_sample_info.csv`


**Pseudobulked snRNA-seq**

* cell type fine: `processed-data/08_pseudoBulkDGE_sn/01_pseudobulk_data_sn/sce_pseudo_DGE-cell_type_anno.RDS`

* cell type broad: `processed-data/08_pseudoBulkDGE_sn/01_pseudobulk_data_sn/sce_pseudo_DGE-cell_type_broad.RDS`


### SRT/Visium

**SpatialExperiment**

* project path: `"processed-data/spe_objects/spe_ERC_annotated"`

* pre-QC version: `processed-data/spe_objects/spe_raw.rds` (Not on github due to large file size: 4G)

```
## Load HD5F spe
spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
# class: SpatialExperiment 
# dim: 30494 122202 

```

** Input files **

* FASTQ: `raw-data/FASTQ`

* Images: `raw-data/FASTQ/Images`

* SAMPLE ID table: `processed-data/02_build_spe/sample_info.csv`


**Pseudobulked Visium data**

`processed-data/09_pseudoBulkDGE_Visium/01_pseudobulk_data_Visium/spe_pseudo_DGE.RDS`


## Internal

JHPCE location: `/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/`

snRNA-seq main data JHPCE path: `"/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/sce_objects/sce_ERC_subcluster"`

SRT main data JHPCE path: `/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/spe_objects/spe_ERC_annotated`



<script type='text/javascript' id='clustrmaps' src='//cdn.clustrmaps.com/map_v2.js?cl=ffffff&w=288&t=n&d=Bh6nY5iydLVTKEEP1j9OciyA69GplP0R5EC0Or6nFLU'></script>
