# Br5529 donor-impact analysis: coherent plain-language narrative

## Core question
Why does removing Br5529 change the APOE carrier analysis from no DEGs to many DEGs?

## Step 0: Observation (not explanation)
- 21 donors: **0 DEGs** (best FDR = **0.1907**)
- Drop Br5529 (20 donors): **138 DEGs** (best FDR = **0.00425**)

## Step 1: Is the focus donor globally unusual in expression space?
- PCA distance rank: **2/21**
- Low-correlation rank: **2/21**
- Extreme-gene burden rank: **5/21**

## Step 2: Could this mainly be a simple covariate artifact?
- Rank convention reminder: **1/N means most outlying** on that metric.
- Mahalanobis rank: **12/21** (not outlying)
- Leverage rank: **5/21** (near outlier)
- chrM ratio rank: **2/21** (strong outlier)

## Step 3: Directional shift on emergent DEGs
- Fraction with positive contrast gain after dropping Br5529: **0.949**
- Median contrast gain: **0.206**
- Focus donor DEG20 median z (within carrier): **1.866**

## Step 4: Same-carrier context
- Focus donor median-z rank on DEG20 genes within E2_carrier donors: **1/10**

## Integrated conclusion
Br5529 shows a coherent expression-space outlier pattern and a donor-direction shift consistent with attenuation of the carrier contrast.
