# Data File paths and code

## 1. Locus Coeruleus (LC) Data (Provided by Bernie)
*Relevant for the "Orange" and "Blue" workflows (Cross-region mediation).*

### **A. LC & Astrocyte Domain DEGs (APOE E4 vs E2)**
**File Path (JHPCE):**
`/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/12_DEanalyses_removedsampsAndFinalNMseg/02b-voomLmFit_objs_and_DE_tables_LabSex_LabE2vE4_LabSingleGenos_LabAnces_LabAncesE4E2_noMixClusts_YRIcontAndCat_JUL2025.RDS`

**Description & Comments:**
*   **Leo's Note:** The original message contained a typo where a `.txt` filename was pasted inside this path. The path above is the corrected version.
*   **Bernie's Description:** This is a multi-layered list containing differential expression (DE) analysis inputs, intermediates, and result tables.
*   **How to access specific results:**
    *   The tables you need are located within the object at: `thisRds$cont_ances$E2E4label$detabs[c("LC_E4vE2","Astro_E4vE2")]`.
    *   These results control for genomic ancestry as a continuous measure.

### **B. Neuromelanin (NM) Intensity DEGs**
**File Path (JHPCE):**
`/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/12_DEanalyses_removedsampsAndFinalNMseg/11-LCNMspots_NMmetric_LMMed_by_geneLogcounts_3apoeFactorizations.RDS`

**Description & Comments:**
*   **Bernie:** This analysis looks at NM pixel intensity controlling for APOE (E4 vs E2 carrier status).
*   **How to access specific results:**
    *   This is a nested list. Relevant tabular results are at: `dat$E4carrierYN$Intensity_NM`.
    *   This table contains coefficients, test stats, and raw Linear Mixed Model (LMM) P-values.
    *   Filter for genes using `dat<-dat[gene_name==variable]`.
    *   **Significance:** P-values are in columns `fdr_gene` (FDR corrected) and `bonf_gene` (Bonferroni corrected).
    *   **Directionality:** A **negative** coefficient value means higher expression is associated with **darker** (stronger) neuromelanin. A positive value is associated with lighter (weaker) neuromelanin.

### **C. NM+ vs NM- Spot DEGs**
**File Path (JHPCE):**
`/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/12_DEanalyses_removedsampsAndFinalNMseg/04-LCNMpos_vs_LCNMneg_DE.txt`

**Description:**
*   Text file containing differential expression results for Neuromelanin positive (NM+) versus Neuromelanin negative (NM-) LC spots.

---

## 2. Entorhinal Cortex (ERC) Data (Provided by Louise)
*Relevant for the "Green" workflow (Within-region mediation).*

### **A. Oligodendrocyte (Oligo.3) DEGs (Carrier Model)**
**File Paths (Relative to project root `LFF_spatial_ERC`):**
*   **RDS:** `processed-data/13_compile_DGE/01_compile_DGE/sn_fine/DGE_results_carrier_sn_fine.Rds`
*   **CSV:** `processed-data/13_compile_DGE/01_compile_DGE/sn_fine/vlmf_model_summary_sn_fine.csv`

**Description & Comments:**
*   **Louise:** These files contain the fine-resolution single-nucleus differential gene expression results.
*   **Usage:** Use the `vlmf` stats. This specific file contains the "Combined all sample 'carrier' model".

### **B. Oligodendrocyte (Oligo.3) DEGs (Ancestry Specific)**
**File Path (Relative to project root `LFF_spatial_ERC`):**
*   **CSV:** `processed-data/13_compile_DGE/05_compile_DGE_ancestry/sn_fine/DGE_results_ancestry_sn_fine.csv`

**Description & Comments:**
*   **Louise's Description:** Contains ancestry-specific model results. European Ancestry (EA) and African Ancestry (AA) stats are noted in the contrast column.

---

## 3. Code & Scripts (Provided by Louise & Leo)

### **A. Global Differential Expression Script (Primary Template)**
**GitHub Path:**
`https://github.com/LieberInstitute/LFF_spatial_ERC/blob/devel/code/12_voomLmFit/01_Clusterwise_voomLmFit.R`

**Comments:**
*   **Louise's Note:** Use this script to find Oligo.3 DEGs by APOE carrier. Note that they only examined the 'carrier' model in this script.
*   **Leo's Instruction:** This is the script you must adapt for the mediation analysis (specifically for the "Green" workflow),.

### **B. Ancestry-Specific Differential Expression Script**
**GitHub Path:**
`https://github.com/LieberInstitute/LFF_spatial_ERC/blob/devel/code/12_voomLmFit/04_Clusterwise_voomLmFit_ancestry.R`

**Comments:**
*   Reference for finding Oligo.3 DEGs specific to each ancestry.

### **C. Data Compilation Script**
**GitHub Path:**
`https://github.com/LieberInstitute/LFF_spatial_ERC/blob/devel/code/13_compile_DGE/01_compile_DGE.R`

**Comments:**
*   Helpful context for how the DGE lists were compiled.

### **D. LC Spatial Domain Data (Integration Reference)**
**GitHub Path:**
`https://github.com/LieberInstitute/MFA_LC-ERC_pilot/blob/81fb5c1d51c3f9a00cba701d6ba0f5975f8a7087/code/00_project_prep/00_data_check.R#L17`

**Comments:**
*   **Leo's note:** This path points to the LC SpD data object.
*   **Louise's note:** This is the object Bernie recommended for the LC-ERC integration, though she notes it is from September 2024 and may need updating.
