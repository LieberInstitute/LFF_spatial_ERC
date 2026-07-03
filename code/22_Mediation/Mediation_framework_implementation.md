# Mediation Framework Implementation

## Purpose

`01_mediation_screening.R` implements an APOE mediation-screening workflow for testing whether candidate mediator gene-expression features attenuate the APOE carrier association observed in ERC Oligo.3 outcome genes. The current script is a screening engine: it refits outcome models with one candidate mediator at a time, records the APOE and mediator coefficients, and compares each joint model against a baseline APOE-only outcome model.

The script does not change the upstream differential-expression analyses that define the outcome genes or mediator candidates. It consumes those upstream results and standardizes them into a common mediator contract.

## Execution and Model Configuration

The script uses the following R objects to configure the analysis:

- `opt`: parsed command-line options from `getopt()`. The most important field is `opt$mediation`, which selects the mediator loader branch.
- `FDRthr`: adjusted P value threshold used when comparing baseline and mediator-adjusted outcome model results.
- `medSelFDR`: adjusted P value threshold used by mediator loaders when selecting candidate mediator genes.
- `drop_donors_global`: global donor exclusions, currently `c("Br1289")`.
- `drop_donors_lc`: LC-specific donor exclusions, set inside LC mediator branches when needed.
- `workers`: number of `BiocParallel::MulticoreParam()` workers.
- `test_rows`: optional row indices parsed from `MEDIATION_TEST_ROWS` for smoke-test runs.
- `forceSingleBaseline`: controls whether all mediator runs reuse the largest donor-set baseline; currently `TRUE`.

The model formulas are also explicit code objects:

```r
baseline_formula_str <- "0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio"
mediation_formula_str <- paste(baseline_formula_str, "+ med_vec")
baseline_formula <- as.formula(paste("~", baseline_formula_str))
mediation_formula <- as.formula(paste("~", mediation_formula_str))
```

`med_vec` is the one-mediator numeric expression vector inserted into the sample metadata by `build_design()` before `model.matrix()` is called.

## Supported Options

The command-line option `--mediation` controls the mediator loader and preserves four current values:

- `erc_astro`: Green design. ERC astrocyte APOE carrier DEGs are used as within-region mediator candidates.
- `lc_astro`: Blue design. LC astrocyte-domain APOE DEGs are used as cross-region mediator candidates.
- `lc_nm`: Orange design. LC NM-positive versus NM-negative associated genes are used as cross-region mediator candidates, then intersected with LC APOE DEGs for the mediator-selection gate.
- `lc_nm_int`: Orange design. LC NM-intensity associated genes are used as cross-region mediator candidates, then intersected with LC APOE DEGs for the mediator-selection gate.

The LC inputs are expected under `processed-data/22_Mediation`. If they are not present locally, the current parent project copy on `gxlin:~/work/R/LFF_spatial_ERC/` has the needed existing input files.

## Core Data Structures

### `sce_pb`

`sce_pb` is the ERC pseudobulk `SingleCellExperiment` loaded from `processed-data/08_pseudoBulkDGE_sn/01_pseudobulk_data_sn/sce_pseudo_DGE-cell_type_anno.RDS`. It supplies:

- the `sce_base` subset for `out_clus <- "Oligo.3"`;
- ERC astrocyte expression for the `erc_astro` mediator option;
- donor metadata such as `BrNum`, `APOE_syn`, `Sex`, `Age`, `Anc_Afr`, `pseudo_expr_chrM_ratio`, and `exp_round`.

`Br1289` is dropped globally before the mediation workflow begins.

### `med_degs`

`med_degs` is the standardized mediator-candidate table. Every mediator option must create a data.table with at least:

- `cluster`: mediator source label, such as `Astro.1`, `LCAstro`, `LC_NMpos_vs_NMneg`, or an analogous new option label;
- `gene_id`: gene identifier used to match expression rows;
- `gene_name`: display gene symbol/name;
- `med_gene_id`: unique mediator row key, conventionally `cluster|gene_id`.

The script then adds:

- `run`: display/run label, `cluster|gene_id|gene_name`;
- `iter`: stable integer identifier used in output filenames.

### `med_expr`

`med_expr` is the mediator expression matrix. Its rows must match `med_degs$med_gene_id`, and its columns must be donor `BrNum` values. Values should be numeric expression measures suitable to enter the ERC outcome model as `med_vec`.

For `erc_astro`, `med_expr` is assembled from ERC astrocyte `logcounts` in `sce_pb`. For `lc_astro`, it is LC astrocyte log-CPM computed from a donor-level `DGEList`. For `lc_nm` and `lc_nm_int`, it is loaded from precomputed LC NM expression matrices.

### `dge_base`

`dge_base` is the ERC Oligo.3 outcome object used for model fitting. After mediator loading, it is restricted to donors with mediator expression columns, filtered with `edgeR::filterByExpr()` using the baseline design, and normalized with `edgeR::calcNormFactors()`.

### Donor Plans

`mediator_plan_dt` records the usable donor set for each mediator row after dropping donors with missing mediator values. Each mediator must have at least 10 donors to be eligible. The current script uses `forceSingleBaseline <- TRUE`, which means the largest available donor set baseline is reused for comparisons when mediator donor sets differ.

### Run Hashes and Caching

The script builds a `run_hash` from:

- mediation option and thresholds;
- model formulas and donor-drop settings;
- input file fingerprints;
- `med_degs`, `med_expr`, and outcome gene hashes;
- outcome donor set.

Each baseline and mediator fit writes a sidecar text file with cache keys. Cached outputs are reused only when these keys match the current run signature.

## Model-Fitting Functions

### `run_baseline()`

For each planned donor set, `run_baseline()` subsets `dge_base` to `common_donors`, builds `des` with `build_design(dge_sub$samples, baseline_formula, batch)`, and fits the model with `voomLmFit()`. `batch` is `exp_round` and is passed as the blocking variable. The APOE carrier coefficient is extracted with:

```r
carrier = "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)"
```

The baseline output table contains the APOE carrier coefficient for all retained Oligo.3 outcome genes. Rows with `adj.P.Val < FDRthr` are stored as the baseline significant gene set for the donor key.

### Mediator Loader Branches

Mediator selection happens before the shared model-fitting code. Each loader branch must create `med_degs` and `med_expr`:

- `erc_astro` selects ERC astrocyte carrier DEGs from compiled ERC DGE results.
- `lc_astro` selects LC astrocyte APOE DEGs from LC spatial-domain DGE results.
- `lc_nm` and `lc_nm_int` start from NM-associated candidate tables and then retain only genes that also pass LC APOE DGE filtering in `LC_E4vE2_DGEout_Bernie.rds`.

That last detail is important: the orange design is biologically motivated by NM association, but the implementation currently enforces APOE association by intersection with the LC APOE DGE table.

### `run_one()`

`run_one()` fits one model per mediator. It resolves one row from `med_degs`, retrieves the planned donor set from `mediator_plan_run`, subsets `dge_base`, and pulls the mediator expression vector from:

```r
med_vec <- as.numeric(med_expr[med_gene_id, donors_i])
```

The design matrix is created from:

```r
design_obj <- build_design(dge_sub$samples, mediation_formula, batch, med_vec = med_vec)
```

The same fitted model is queried twice with `topTable()`:

- once for the APOE carrier contrast;
- once for the mediator coefficient `med_vec`.

The summary logic marks outcome genes as attenuated when they were significant in the baseline APOE model but are no longer significant for APOE after including `med_vec`. A mediator/outcome pair passes both practical screening criteria when the outcome is attenuated and the mediator coefficient is significant in the joint model.

## Outputs

Important output files include:

- `mediation_mediator_index.tsv`: standardized mediator index for the latest run.
- `baseline_index.tsv`: baseline model outputs by donor set.
- `mediation_run_index.tsv`: one row per completed mediator model with output paths and hashes.
- `mediation_vs_baseline_summary.tsv.gz`: attenuation summary by mediator.
- `run_history.tsv`: append-only run metadata.
- `*.tab.gz`: APOE carrier coefficient tables for individual mediator models.
- `*.M.tab.gz`: mediator coefficient tables for individual mediator models.
- `*.loss.txt`, `*.gain.txt`, `*.mediated.lst`: optional gene-list side outputs when non-empty.

## Adding a New Mediator Dataset Option

A new mediator option should be added only by inserting a new loader branch that produces the existing contract. The core model-fitting functions should not need to change if the new dataset can provide donor-level mediator expression.

Checklist:

1. Add the option name to the help text, validation list, and `fdr_defaults`.
2. Add a loader branch before the final "not implemented" stop.
3. Define the mediator-selection rule, such as APOE-associated DEGs, NM-associated genes intersected with APOE DEGs, or another documented upstream screen.
4. Load or build `med_degs` with `cluster`, `gene_id`, `gene_name`, and `med_gene_id`.
5. Load or build `med_expr` as a numeric matrix with rownames matching `med_degs$med_gene_id` and columns named by donor `BrNum`.
6. Restrict columns to donors shared with ERC Oligo.3 outcome data, applying any option-specific donor exclusions.
7. Confirm that each mediator row has enough non-missing donor values for the minimum donor threshold.
8. Add all input files to `run_input_files` so `run_hash` and cache invalidation remain valid.
9. Use `MEDIATION_TEST_ROWS` for a small smoke test before launching a full run.
10. Document the upstream source model and any contrast differences from the ERC outcome model.

Minimum mediator contract:

```text
med_degs:
  cluster
  gene_id
  gene_name
  med_gene_id

med_expr:
  numeric matrix
  rownames == med_degs$med_gene_id
  colnames == donor BrNum values
```

The most common failure modes are mismatched gene identifiers, duplicated `med_gene_id` values, missing donor columns, and mediator rows that pass the upstream DEG filter but are removed by expression filtering.

## Supporting a Different Outcome Dataset

Changing the base mediated/outcome gene dataset is a larger change than adding a mediator option. It would affect:

- `out_clus`, currently fixed to `Oligo.3`;
- the source object used for `sce_base` and `dge_base`;
- donor metadata columns and covariates;
- the baseline model formula and contrast;
- the baseline DEG definition used for attenuation;
- output labeling and gene-map assumptions.

If the new outcome still lives in the same ERC pseudobulk SCE and uses the same APOE carrier contrast, this may be a modest parameterization change. If it comes from a different region, platform, cell type system, or genotype encoding, the outcome setup should be refactored into an explicit outcome loader parallel to the mediator loaders.

## Xenium/LR Follow-Up

Louise's `code/21_Xenium/18_xenium_DEG_mediation_LR.R` script is best treated as downstream validation/context for mediation-derived mediator/outcome pairs. It checks whether selected mediator and outcome genes are represented in Xenium DGE and ligand-receptor resources. It is not part of the core mediation refit loop and should not be described as defining mediation hits.
