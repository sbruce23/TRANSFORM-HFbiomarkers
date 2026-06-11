# TRANSFORM-HF Biomarker Study — Analysis Code Repository

This repository contains annotated code supporting the analyses, figures, and tables in the manuscript:

> **Biomarker-Defined Heterogeneity in NT-proBNP Response to Torsemide Versus Furosemide After Heart Failure Hospitalization: Results from the TRANSFORM-HF Ancillary Biomarker Study**
> Lauren Cooper, MD · Scott A. Bruce, PhD · Stephen Seliger, MS, MD · *et al.*

The manuscript evaluates whether baseline biomarker profiles identify differential biological responses to loop diuretic strategy (torsemide vs. furosemide) in patients discharged after a heart failure hospitalization.

> ⚠️ **Note:** Individual-level participant data are not included in this repository. The code is provided to document the analytic workflow and support adaptation of these methods in future work.

---

## Study Overview

The TRANSFORM-HF Ancillary Biomarker Study enrolled **305 participants** (163 from the TRANSFORM-HF randomized trial and 142 from an observational extension cohort) discharged on either torsemide or furosemide following acute heart failure hospitalization. Biomarkers were measured at enrollment and at 90 days post-discharge.

**Key design challenges addressed:**

1. **Non-random treatment assignment** in the observational extension cohort — addressed using Inverse Probability of Treatment Weighting (IPTW).
2. **Differential dropout** (~35% of participants did not return for 90-day blood collection) — addressed using Inverse Probability of Censoring Weighting (IPCW).

A combined inverse probability weight (IPW) — the product of the treatment and dropout weights — was applied in all regression analyses.

---

## Repository Contents

### 1. Main Analysis: Tables and Figures

📄 **[`Transform_Analysis_GH.html`](https://sbruce23.github.io/TRANSFORM-HFbiomarkers/Transform_Analysis_GH.html)**

This rendered R Markdown document reproduces all tables and figures from the manuscript, including:

- **Baseline characteristics** (Table 1, unweighted; Supplemental Table 2, IPW-weighted)
- **IPW-weighted linear regression models** for 90-day change in NT-proBNP (log₂ scale) — with biomarker-by-treatment interaction terms for NT-proBNP, GDF-15, hs-CRP, hs-cTnT, hs-cTnI, Cystatin C, and Creatinine (Table 2, Figure 2)
- **IPW-weighted logistic regression models** for the probability of achieving NT-proBNP < 1,000 pg/mL at 90 days — with estimated odds ratios at the 25th and 75th percentile of each baseline biomarker (Table 3, Figure 3)
- **Bootstrap confidence intervals** for predicted outcomes at the 25th and 75th percentiles of each baseline biomarker
- **Causal tree (CART) analyses** identifying biomarker-defined subgroups with heterogeneous diuretic treatment effects (Figure 4)

The analysis is implemented in **R**, using IPW-weighted regression, bootstrap resampling, and causal tree methods.

---

### 2. IPW Model Notebooks

The inverse probability weights used in the main analysis were generated using **XGBoost** classification models implemented in **Python**. Two separate propensity score models were fit:

---

#### 2a. Dropout (IPCW) Model

📓 **[`explanatory_dropout_model_GH.ipynb`](https://github.com/sbruce23/TRANSFORM-HFbiomarkers/blob/main/explanatory_dropout_model_GH.ipynb)**

This notebook estimates each participant's **probability of dropping out** (i.e., not providing a 90-day blood sample), which is used to construct Inverse Probability of Censoring Weights (IPCW).

**Model details:**

| Item | Description |
|---|---|
| **Outcome** | Binary indicator of dropout (1 = no 90-day follow-up) |
| **Training data** | 67% of participants (n = 191); held-out 33% (n = 90) |
| **Algorithm** | XGBoost classifier (`xgboost` library) |
| **Class imbalance** | Minority class (dropout) upweighted via `scale_pos_weight` |
| **Hyperparameter tuning** | Grid search with predefined train/holdout split; optimized F1 score |
| **Hyperparameter grid** | `n_estimators`, `colsample_bytree`, `max_depth`, `subsample`, `min_child_weight` |
| **Covariates** | Site, demographics (age, sex, race, ethnicity), cardiovascular history (diabetes, CAD, CKD, hypertension, prior HF, smoking, etc.), insurance type, baseline biomarkers (NT-proBNP, GDF-15, hs-CRP, Troponin I, TNT, Cystatin C, Creatinine, COV2-S, COV2-N), treatment group, and randomization indicator |

**Notebook outputs include:** predicted dropout probabilities, KDE plots of predicted probabilities by observed dropout status, AUC training curve, classification metrics (accuracy, sensitivity, specificity, F1, AUC, Brier score), and feature importance plot.

---

#### 2b. Treatment (IPTW) Model

📓 **[`explanatory_treatment_model_GH.ipynb`](https://github.com/sbruce23/TRANSFORM-HFbiomarkers/blob/main/explanatory_treatment_model_GH.ipynb)**

This notebook estimates each **observational-cohort participant's propensity score** — the conditional probability of receiving torsemide versus furosemide given baseline covariates — which is used to construct Inverse Probability of Treatment Weights (IPTW).

> **Note:** Participants from the **randomized** component of the study were assigned a fixed treatment probability of 0.5, reflecting the 1:1 randomization scheme. The XGBoost model was applied only to the **observational extension cohort** (n = 142), where diuretic assignment was not randomized.

**Model details:**

| Item | Description |
|---|---|
| **Outcome** | Binary treatment indicator (1 = torsemide, 0 = furosemide) |
| **Training data** | Observational cohort participants only (n = 90 training, n = 45 held-out) |
| **Algorithm** | XGBoost classifier (`xgboost` library) |
| **Class imbalance** | Minority class (torsemide) upweighted via `scale_pos_weight` |
| **Hyperparameter tuning** | Grid search with predefined train/holdout split; optimized F1 score |
| **Hyperparameter grid** | Same candidate grid as the dropout model |
| **Covariates** | Site, demographics, cardiovascular history, insurance type, and baseline biomarkers (same set as dropout model, excluding treatment group and randomization indicator) |

**Notebook outputs include:** predicted treatment probabilities, KDE plots of predicted probabilities by observed treatment, AUC training curve, classification metrics, and feature importance plot.

---

### Combined Weights

For each participant who completed the 90-day follow-up, the final IPW was computed as:

$$w_i = \frac{1}{P(T_i \mid X_i)} \cdot \frac{1}{P(R_i = 1 \mid X_i)}$$

where $T_i$ is the observed diuretic assignment and $R_i = 1$ indicates the participant provided a 90-day follow-up blood sample. These combined weights were applied in all regression analyses in the main R analysis document.

---

## Software & Dependencies

### R (Main Analysis)
- R ≥ 4.3.0
- Key packages:
  - `ggplot2`, `ggpubr`, `gridExtra` — plotting
  - `tidyverse` — data manipulation
  - `boot` — bootstrap confidence intervals
  - `survey` — IPW-weighted summary statistics
  - `Hmisc` — weighted means
  - `broom`, `modelsummary` — model output formatting
  - `gtools` — p-value stars
  - `htetree` — causal tree estimation (`causalTree`, `estimate.causalTree`)
  - `rpart`, `rpart.plot` — tree fitting and visualization
  - `partykit` — party-based tree predictions
  - `glmnet`, `GGally`, `corrplot` — exploratory modeling and correlation analysis

### Python (IPW Models)
- Python ≥ 3.10
- Key packages:
  - `xgboost` — XGBoost classifier
  - `scikit-learn` — `GridSearchCV`, `StratifiedKFold`, `PredefinedSplit`, metrics
  - `pandas`, `numpy` — data handling
  - `matplotlib`, `seaborn` — visualization
  - `tqdm`, `pickle` — progress tracking and model serialization
- GPU training used during model development (`device='cuda'`); CPU execution is supported by removing the `device` argument

---

## Citation

If you use this code, please cite the manuscript (citation details to be added upon publication).

---

## Questions

For questions about the code or methodology, please open a GitHub issue or contact the repository maintainer.
