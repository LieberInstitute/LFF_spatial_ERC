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
For the 3 mediation designs discussed, `X` and `Y` remain the same. The designs differ by mediator pool and donor set.

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
**Workflow: Baron & Kenny Mediation Screening** : being prototyped with **ERC Astrocyte** genes as mediators but is designed to be applied to **Locus Coeruleus (LC)** features (NM-linked signals or LC Astrocyte DEGs) later.

**Variables:**
* **X (Predictor):** APOE Genotype (E4+ vs. E2+).
* **Y (Outcome):** ERC Oligo.3 Gene Expression/DEGs (1,022 Baseline DEGs).
* **M (Mediator):** Candidate Gene Expression/DEGs (e.g., from LC or ERC Astrocytes).

#### **Step 1: Establish the Total Effect (X -> Y)**

* **Goal:** Confirm that APOE genotype significantly affects the outcome (Oligo.3 expression) *before* accounting for any mediator.
* **Action:** Run the baseline differential expression model: Y ~ X + Covariates
* **Criteria:**
    - The Oligo.3 gene must be differentially expressed (DEG).
    - **Threshold:** APOE coefficient **FDR < 0.05**.
    - *Note: The current list of 1,022 Oligo.3 DEGs already satisfies this step.*

#### **Step 2: Establish the Mediator's Relationship to X (X -> M)**

* **Goal:** Confirm that APOE genotype significantly affects the candidate mediator gene.
* **Action:** Select only mediator candidates that are themselves DEGs in their source cell type (e.g., Astro.[1-3]): M ~ X + Covariates
* **Criteria:**
    - The candidate mediator gene must be a significant DEG.
    - **Threshold:** APOE coefficient **FDR < 0.05**.
    - *Note: This was enacted by selecting the Astrocyte or LC DEGs as mediation candidates.*

#### **Step 3: Test Mediation (Joint Model: X + M -> Y)**

* **Goal:** Determine if M affects Y *and* if it explains the original effect of X on Y.
* **Action:** Fit the joint model with the mediator expression vector included as a covariate:
      Y ~ X + M + Covariates
* *Implementation:* Besides the APOE contrast, defined the M contrast as well. Then run `topTable()` twice on this single model fit: once for the Mediator coefficient and once for the APOE coefficient.

* **Criteria (Must pass both):**
    - **Step 3a ( M -> Y significance):** The mediator must significantly predict the outcome.
        - **Threshold:** Mediator coefficient **FDR < 0.05**.
    - **Step 3b ( X -> Y attenuation):** The effect of APOE on the outcome should lose significance (indicating the effect was "mediated" by M).
        - **Threshold:** APOE coefficient **FDR >= 0.05** (i.e., statistically indistinguishable from 0).

### Practical (R) notes

The base design formula for the Oligo.3 cell type differential gene analysis (Y~ X + covariates) on the pseudobulk SingleCellExperiment (ScE) object `spe_pb`:

 ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio

`APOE_syn` is a factor with these levels and distributions in the ScE colData():

```
> table(colData(sce_pb)$APOE_syn)

E2.E2 E2.E3 E3.E4 E4.E4
  162   252   289   180
```

The APOE contrast is defined like this:
```
 makeContrasts(carrier = "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)", ...)
```
We add M as a covariate to this initial Y ~ X design formula by extracting a `med_vec` numeric vector with the gene expression for the candidate mediator genes (which are significant DEGs in other DGE analyses), so the formula for Step 3. becomes:

~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio + med_vec
