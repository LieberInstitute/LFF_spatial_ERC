## LC <-> ERC mediation-screening framework (APOE risk mechanisms)

In this project we are performing a mediation analysis on potential mechanisms of Alzheimer disease, investigating how changes shown in a Locus Coeruleus (LC) study (https://www.biorxiv.org/content/10.1101/2025.10.29.685354v1.full), could impact (mediate) changes in Entorhinal Cortex (ERC) as presented in another study (https://www.biorxiv.org/content/10.1101/2025.11.20.689483v1.full).

This is an practical implementation of a "mediation analysis" framework to connect **APOE genotype (E4 vs E2 carrier status)** to **ERC Oligo.3 expression phenotypes** (outcome), through potential mediation from **Locus Coeruleus (LC) gene expression** (including neuromelanin NM-linked signals).

A key practical enabler is the 24 donors overlap between the ERC and LC studies supporting this cross-region mediation modeling.

### Key idea (operational "mediation" screening):
* **X** = APOE carrier status (E4+ vs E2+)
* **Y** = ERC Oligo.3 gene expression outcomes (often starting from ~1,000 APOE-associated Oligo.3 DEGs).
* Testing candidate **M** (mediators) *one at a time* by re-fitting outcome models:
    * Baseline: `Y ~ X + covariates`
    * Adjusted: `Y ~ X + M + covariates`
* A mediator "hit" is when including **M** in the regression **attenuates** the APOE effect on **Y** (e.g., beta coefficient shrinks toward 0 / loses significance) across a meaningful subset of Oligo.3 outcomes.

### The 3 mediation designs discussed in the video (differ by mediator pool and donor set)
(these are color-coded in the video presentation)
**1. Design 1 (Orange):**
**Cross-region LC domain / NM-driven mediators => ERC Oligo.3 outcomes (n = 24 overlap donors)**
* **X:** APOE (E4+ vs E2+)
* **M:** LC spatial-domain gene features, including:
    * APOE-associated LC domain DEGs
    * NM+ vs NM- DEGs
    * NM-intensity association genes
* **Y:** ERC Oligo.3 outcomes (gene-by-gene over the Oligo.3 APOE DEG set)
* **Goal:** test whether LC state/NM-linked biology can explain part of the APOE signal seen in ERC oligodendrocytes.

**2. Design 2 (Blue):**
**Cross-region LC astrocyte-domain mediators => ERC Oligo.3 outcomes (n = 24 overlap donors)**
* **X:** APOE (E4+ vs E2+)
* **M:** APOE-associated DEGs from the **LC astrocyte spatial domain** (astro/proximity signal near LC neurons?) (use their expression as mediators)
* **Y:** ERC Oligo.3 outcomes
* **Goal:** a more focused LC=>ERC hypothesis: APOE-driven astro changes in LC relate to (or mediate) APOE-driven Oligo.3 changes in ERC.

**3. Design 3 (Green): ERC-only astrocyte mediators => ERC Oligo.3 outcomes (n = 30 ERC donors)**
* **X:** APOE (E4+ vs E2+)
* **M:** APOE-associated DEGs from ERC astrocyte clusters (e.g. ~88 in the initial plan)
* **Y:** ERC Oligo.3 outcomes
* **this is a within-region mediation screen (no LC involved)**
    * initial prototyping run/tooling
    * testing if APOE expression in ERC astrocytes mediate/affect ERC Oligo.3 APOE expression

### Practical implementation plan (high level):
* Start with the "green" design #3: **main Oligo.3 APOE DEG list** ; later, consider ancestry-stratified or APOE x ancestry interaction analyses
* run one mediator gene at a time (array job)
* reuse the existing ERC DE pipeline (limma->voomLmFit modeling), only adding the mediator gene expression as an extra covariate in the model, one at a time
* aggregate the many topTable outputs into a **Oligo.3 outcome genes x mediator genes** summary:
    * how often APOE remains significant after adjusting for M
    * effect attenuation metrics (e.g., beta changes for APOE)
