# GLUNEX — Multimodal Type 2 Diabetes Risk Prediction System

> Integrating Clinical, Lifestyle, and Gene Expression Data via Late Fusion

**Author:** Talha Rehman  
**Degree:** Master of Science  
**Submitted:** March 2026

---

## Overview

**Glunex** is a multimodal machine learning system for Type 2 Diabetes (T2D) risk prediction. It integrates three independent data modalities — structured clinical data with free-text clinical notes, self-reported lifestyle survey data, and blood transcriptomic gene expression data — through a principled late-fusion architecture.

T2D affects approximately 589 million adults globally (2024) and is projected to reach 853 million by 2050. Existing prediction models predominantly rely on single-domain clinical data, failing to capture the full biological complexity of T2D. Glunex addresses this by fusing complementary signals from clinical, behavioural, and molecular dimensions.

---

## System Architecture

Glunex uses a **late-fusion architecture**: three base models are trained independently on separate datasets, and their calibrated output probabilities are combined through a logistic regression meta-model.

```
Clinical Data (EHR + Notes)   →   Clinical Model (XGBoost)     →  p_clinical
Lifestyle Data (CDC BRFSS)    →   Lifestyle Model (LightGBM)   →  p_lifestyle   →  Fusion Meta-Model  →  Final Risk Score
Gene Expression (Blood RNA)   →   Gene Expression Model (XGB)  →  p_gene
```

The system supports **four input scenarios** depending on data availability:
1. All three modalities (clinical + lifestyle + gene expression)
2. Clinical + lifestyle only
3. Clinical only
4. Lifestyle only

---

## Model Performance Summary

| Modality | Algorithm | Dataset | AUC |
|---|---|---|---|
| Clinical | XGBoost | Kaggle Diabetes Dataset (99,986 patients) | 0.9754 |
| Lifestyle | LightGBM | CDC BRFSS Survey (253,680 respondents) | 0.8241 |
| Gene Expression (Blood) | XGBoost | GSE184050 (116 samples, LOO-CV) | 0.7679 |
| **Fusion Meta-Model** | **Logistic Regression** | **Synthetic (10,000 samples)** | **0.9897** |

**Fusion Brier Score:** 0.0272  
**SHAP Modality Attribution:** Clinical 62.7% · Lifestyle 22.4% · Gene Expression 14.9%

---

## Repository Structure

```
Multimodal-T2D-Classifier-Model-Training/
│
├── README.md
│
├── clinical_model/
│   └── model_train.ipynb          # Full pipeline: preprocessing, NLP, XGBoost training, SHAP, calibration
│
├── lifestyle_model/
│   └── lifestyle_model.ipynb      # Full pipeline: BRFSS preprocessing, feature engineering, LightGBM, calibration
│
├── fusion_model/
│   └── late_fusion.ipynb          # Percentile-rank normalisation, meta-model training, SHAP attribution
│
└── gene_expression_model/
    ├── blood_model/               # R scripts for blood transcriptomics (GSE184050)
    │   └── *.R                    # Step-by-step: data loading → QC → WGCNA → classification → validation
    └── islet_model/               # R scripts for islet transcriptomics (GSE50244)
        └── *.R                    # Step-by-step: data loading → QC → WGCNA → classification → validation
```

---

## Models In Detail

### 1. Clinical Model (`clinical_model/`)

- **Dataset:** Synthetic EHR dataset based on the Kaggle Diabetes Prediction Dataset (99,986 patients, 8.5% diabetic)
- **Features:** Demographics, BMI, HbA1c, blood glucose, hypertension, heart disease, smoking history + free-text clinical notes
- **Feature Engineering:** `age_group`, `bmi_category`, `comorbidity_count`, `metabolic_index`
- **NLP Pipeline:** TF-IDF vectorisation (167 terms) → Truncated SVD (30 components, 97.74% variance explained)
- **Final Feature Matrix:** 54-dimensional (24 structured + 30 text features)
- **Class Imbalance:** Handled with SMOTE on training set only
- **Best Model:** XGBoost (Val AUC = 0.9785, F1 = 0.8060)
- **Calibration:** Isotonic regression; calibrated mean probability (0.0842) closely matches class prevalence (8.5%)
- **Explainability:** SHAP TreeExplainer; top features: HbA1c, metabolic_index, blood glucose

---

### 2. Lifestyle Model (`lifestyle_model/`)

- **Dataset:** CDC Behavioral Risk Factor Surveillance System (BRFSS), 253,680 respondents, 22 features
- **Target:** Binary diabetes risk (prediabetes + diabetes = positive class; 15.76% positive rate)
- **Feature Engineering:** 9 composite features including `Health_Behavior_Score`, `Comorbidity_Count`, `Lifestyle_Risk_Index`, `Age_BMI_Interaction`
- **Class Imbalance:** SMOTE within ImbPipeline (training only, zero leakage)
- **Best Model:** LightGBM (Val AUC = 0.8243)
- **Test Performance:** AUC = 0.8241, Brier Score = 0.1060, Optimal threshold = 0.246
- **Calibration:** Isotonic regression; calibrated mean (0.1578) vs actual prevalence (0.1576)

---

### 3. Gene Expression Model (`gene_expression_model/`)

The most methodologically novel component, comprising two complete transcriptomic classifiers.

#### Blood Model — GSE184050
- **Data:** 116 whole blood RNA-seq samples (50 T2D, 66 controls); DESeq2-normalised counts
- **Preprocessing:** Gene-level QC, log2(x+1) transformation, top 5,000 genes by MAD
- **Network:** Signed hybrid WGCNA using biweight midcorrelation (bicor); soft-thresholding power selected at R² ≥ 0.85
- **Classification:** XGBoost on module eigengenes; LOO-CV AUC = 0.7679
- **External Validation:** Same-tissue (GSE9006) AUC = 0.771 · Cross-tissue islet (GSE76894) AUC = 0.635

#### Islet Model — GSE50244
- **Data:** 62 pancreatic islet samples; same preprocessing pipeline
- **Classification:** XGBoost on WGCNA module eigengenes; Training AUC = 0.836
- **External Validation:** Same-tissue (GSE164416) AUC = 0.608 · Cross-tissue blood AUC = 0.587

#### Golden Gene Analysis
- 28 genes identified as hub genes simultaneously in both blood and islet networks
- Provides a cross-tissue T2D biomarker panel with potential translational value

---

### 4. Fusion Meta-Model (`fusion_model/`)

- **Input:** Calibrated probability outputs from all three base models (p_clinical, p_lifestyle, p_gene)
- **Normalisation:** Percentile-rank normalisation applied to each modality's probability before fusion
- **Training Data:** 10,000 synthetic patients constructed by independent sampling from modality-specific probability pools
- **Meta-Learner:** Logistic Regression
- **Performance:** AUC = 0.9897, Brier Score = 0.0272
- **Explainability:** SHAP attribution — Clinical 62.7%, Lifestyle 22.4%, Gene Expression 14.9%

---

## Datasets

The datasets used in this project are **publicly available** and should be downloaded separately:

| Dataset | Source | Description |
|---|---|---|
| Kaggle Diabetes Prediction Dataset | [Kaggle](https://www.kaggle.com/) | Clinical EHR data (99,986 patients) |
| CDC BRFSS Health Indicators | [Kaggle / CDC](https://www.cdc.gov/brfss/) | Lifestyle survey data (253,680 respondents) |
| GSE184050 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184050) | Blood RNA-seq training data (116 samples) |
| GSE9006 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9006) | Blood RNA-seq external validation |
| GSE50244 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50244) | Islet RNA-seq training data (62 samples) |
| GSE164416 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164416) | Islet RNA-seq external validation |
| GSE76894 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76894) | Cross-tissue islet validation |

> Raw data files are not included in this repository due to size constraints.

---

## Dependencies

### Python (Notebooks)
Key libraries used across the clinical, lifestyle, and fusion notebooks:
```
xgboost
lightgbm
scikit-learn
imbalanced-learn
shap
pandas
numpy
nltk
scipy
matplotlib
seaborn
joblib
```

Install via:
```bash
pip install xgboost lightgbm scikit-learn imbalanced-learn shap pandas numpy nltk scipy matplotlib seaborn joblib
```

### R (Gene Expression Scripts)
Key R packages used:
```
WGCNA
DESeq2
GEOquery
xgboost
caret
ggplot2
dplyr
```

Install via:
```r
install.packages(c("WGCNA", "xgboost", "caret", "ggplot2", "dplyr"))
BiocManager::install(c("DESeq2", "GEOquery"))
```

---

## How to Run

### Clinical Model
Open `clinical_model/model_train.ipynb` and run all cells sequentially. The notebook covers the full pipeline from raw data loading to calibrated probability output and SHAP analysis.

### Lifestyle Model
Open `lifestyle_model/lifestyle_model.ipynb` and run all cells sequentially. Full pipeline from BRFSS data loading to calibrated LightGBM output.

### Gene Expression Model
Run the R scripts in `gene_expression_model/blood_model/` and `gene_expression_model/islet_model/` **in numbered order** (e.g., `01_`, `02_`, etc.). Each script corresponds to a distinct stage of the pipeline.

### Fusion Model
Open `fusion_model/late_fusion.ipynb` after base model outputs are generated. This notebook combines the three calibrated probabilities and trains the meta-learner.

---

## Key Findings

- Multimodal fusion (AUC = **0.9897**) substantially outperforms all individual modality models
- Blood transcriptomic signatures generalise to independent same-tissue cohorts (AUC = 0.771) and retain partial signal in cross-tissue islet validation (AUC = 0.635)
- **28 Golden Genes** were identified as hub genes in both blood and islet co-expression networks, representing potential cross-tissue T2D biomarkers
- SHAP attribution confirms that clinical data dominates prediction (62.7%) while lifestyle and gene expression contribute meaningful complementary signal

---

## Citation

If you use this work, please cite:

> Rehman, T. (2026). *GLUNEX: Multimodal Type 2 Diabetes Risk Prediction System Integrating Clinical, Lifestyle, and Gene Expression Data via Late Fusion.* MSc Thesis.

---

## License

This repository is for academic and research purposes. Dataset usage is subject to the terms of their respective sources (GEO, CDC, Kaggle).
