## LC <-> ERC mediation-screening framework (APOE risk mechanisms)

We discussed a practical "mediation analysis" workflow to connect APOE genotype (E4 vs E2 carrier status) to ERC Oligo.3 expression phenotypes (DGE outcome), through potential mediation/attenuation from Locus Coeruleus (LC) gene expression (LC spatial domain neuromelanin NM-linked signals, and LC astrocytes domain).

A key practical enabler is the 24 donors overlap between the ERC and LC studies supporting this cross-region mediation modeling.

### Key concepts (operational "mediation" screening)

 * `X` = APOE carrier status (E4+ vs E2+)
 * `Y` = ERC Oligo.3 DGE outcomes (baseline DEGs at FDR<.05 starting from 1022 Oligo.3 DEGs for all 30 subjects).
 * Testing candidate M (mediators) gene expression data, one-at-a-time by re-fitting outcome models:
    - Baseline was: `Y ~ X + covariates`
    - Mediator-adjusted: `Y ~ X + M + covariates` (where `M` is the potential mediator gene expression)
 A mediator "hit" is when including `M` in the regression attenuates the APOE effect on `Y` (e.g., beta coefficient shrinks toward 0, DEGs lose significance), across a subset of Oligo.3 DGE outcomes.

### 3 mediation designs

For the 3 mediation designs discussed in the video, `X` and `Y` remain the same. The designs differ by mediator pool and donor set.

#### 1. Design 1 (Orange):
 Cross-region LC domain / NM-driven mediators => ERC Oligo.3 outcomes (n = 24 overlap donors)
 * M: expression of LC spatial-domain gene features:
   - APOE-associated LC domain DEGs (NM+/NM- DEGs + NM-intensity-associated)
 Goal: test whether LC state/NM-linked biology can explain part of the APOE signal seen in ERC Oligo.3.

#### 2. Design 2 (Blue):
 Cross-region LC astrocyte-domain mediators => ERC Oligo.3 outcomes (n = 24 overlap donors)
 * M: expression of APOE-associated DEGs from the LC astrocyte spatial domain
 Goal: a more focused LC=>ERC hypothesis: APOE-driven astro changes in LC relate to (or mediate) APOE-driven Oligo.3 changes in ERC.

#### 3. Design 3 (Green):
 ERC-only astrocyte mediators => ERC Oligo.3 outcomes (n = 26 donors for Astro.1, 30 donors for Astro.2/Astro.3)
 * M: expression of APOE-associated DEGs from ERC astrocyte clusters (e.g. ~88 in the initial plan).
 Goal: testing if APOE expression in ERC astrocytes mediate/affect ERC Oligo.3 APOE expression (DEGs)
 This is a within-region mediation screen (no LC data involved), also serving as an initial prototyping run/tooling check.

### Practical implementation plan (high level)

* start with the "green" design #3: main Oligo.3 APOE DEG list ; later, consider ancestry-stratified or APOE x ancestry interaction analyses
* baseline carrier design remains the same throughout the 3 scenarios, only donor set changes
* reuse the existing ERC DE pipeline (limma->voomLmFit modeling), only adding the mediator gene expression as an extra covariate in the model, one M gene at a time (array job/BiocParallel loop)
* for each potential mediator (M gene expression covariate ), check the gain/loss in ERC Oligo.3 DEGs outcome.
  - initial screening: DEG loss/gain reporting -- how many baseline DEGs lose significance (FDR<0.05) after adjusting for M gene expression.
