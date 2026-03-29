# Topological Cancer Signature Pipeline (TDA-Based)

This repository provides a complete pipeline for identifying cancer-specific topological features, stable samples, and biologically significant genes from scRNA-seq data.

The workflow integrates R (Seurat-based preprocessing) and Python (Topological Data Analysis + Machine Learning).

---

## STEP 0: REQUIREMENTS

Install the following:

R Packages:

* Seurat
* tidyverse
* data.table
* Hmisc
* Matrix
* readxl
* org.Hs.eg.db
* AnnotationDbi

Python Packages (auto-installed in scripts):

* pandas
* numpy
* matplotlib
* gudhi
* persim
* scipy
* networkx
* openpyxl
* scikit-learn
* mygene

---

## STEP 1: GENERATE CORRELATION MATRIX (R)

Run the R function:
FindCorrMat.R

Function:
FindCorrMat()

Purpose:

* Load scRNA-seq data
* Perform QC, normalization, PCA
* Select top genes
* Build gene correlation network

Output:

* CancerCorrMat.csv

Reference:
See function in code file 

---

## STEP 2: EXTRACT TOPOLOGICAL FEATURES (PYTHON)

Run:
Extract_Topological_Features()

Purpose:

* Construct Vietoris-Rips complex
* Compute persistence diagrams
* Extract TDA features
* Identify highly correlated (significant) genes

Outputs:

* Topological_Features.xlsx
* SignificantGenes.csv
* Betti plots, persistence diagram, landscape plots

Reference:
See function in code file 

---

## STEP 3: COMBINE ALL SAMPLES

Goal:

* Combine outputs from multiple samples

Action:

* Merge all Topological_Features.xlsx into a single file:
  → TDA-Features.xlsx

Each row should contain:

* Sample features
* Cohort label (e.g., Breast, Lung, Healthy Breast, etc.)
* Dataset number

---

## STEP 4: MULTICLASS CLASSIFICATION

Run Python script:
Multiclass classification pipeline

Purpose:

* Classify cancer types using TDA features
* Identify stable samples

Method:

* Polynomial Logistic Regression
* Repeated train-test splits (n=100)

Outputs:

* Precision, Recall, F1-score per cancer
* Stable samples (correctly classified in ALL runs)

Stored variables:

* stable_samples
* stable_samples_by_cancer

---

## STEP 5: BINARY CLASSIFICATION (CANCER vs HEALTHY)

Run Python script:
Binary classification pipeline

User Input:

* cancer_type = "Breast" (or any cancer)

Purpose:

* Compare Stable Cancer vs Healthy samples
* Identify refined stable samples

Method:

* Repeated Stratified K-Fold
* 100% specificity threshold

Outputs:

* Mean accuracy
* Samples contributing to mean accuracy

---

## STEP 6: INTEGRATE TOPOLOGICALLY SIGNIFICANT GENES

Run Python script:
Gene integration pipeline

Purpose:

* Convert Ensembl IDs → Gene Symbols
* Identify genes present in ≥50% of samples

Output:

* Topologically significant cancer-specific gene list (CSV)

Example:

* Breast__Top_Sig_Genes.csv

---

## STEP 7: GENE EXPRESSION PROFILE (R)

Run:
FindGeneProfile.R

Function:
FindGeneProfile()

Purpose:

* Extract expression of topologically significant genes
* Match genes with dataset
* Generate expression profile for DEA

Outputs:

* Recurring genes file
* Gene expression profile

Reference:
See function in code file 

---

## FINAL OUTPUTS

1. CancerCorrMat.csv
2. Topological_Features.xlsx
3. TDA-Features.xlsx
4. Stable samples per cancer
5. Binary classification results
6. Topologically Significant genes (≥50%)
7. Gene expression profiles (for DEA)

---

## PIPELINE SUMMARY

scRNA-seq data
↓
FindCorrMat (R)
↓
CancerCorrMat.csv
↓
Extract_Topological_Features (Python)
↓
Topological Features + Significant Genes
↓
Multiclass Classification
↓
Stable Samples
↓
Binary Classification
↓
Refined Stable Samples
↓
Topologically Significant Gene Integration (≥50%)
↓
FindGeneProfile (R)
↓
DEA-ready Gene Profiles

---

## NOTES

* Always keep consistent naming for cancer types:
  Example: "Breast", "Healthy Breast"

* Ensure:

  * Same gene ID format across datasets
  * Proper file paths

* Stable samples are highly reliable biological signals.

---

## END
