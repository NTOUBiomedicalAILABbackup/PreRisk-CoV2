![PreRisk-CoV2](logo_v2.png)

## Overview

PreRisk-CoV2 is a machine learning framework for pre-exposure risk assessment of SARS-CoV-2 susceptibility using Serum protein biomarkers. The main function is to predict infection risk **before exposure** based on a 7-protein panel identified through K-Nearest Neighbors (KNN) combined with Genetic Algorithm (GA) feature selection. 



If you have any trouble installing or using PreRisk-CoV2, you can post an issue or directly email us. We welcome any suggestions.

---


## Quick Install

*Note*: We suggest you install all packages using conda ([Anaconda](https://anaconda.org/)).

### Prepare the Environment

#### 1. First-time Setup


Create conda environment with required dependencies
```bash
conda create -n PreRisk_CoV2 python=3.9 -y
conda activate PreRisk_CoV2
```


Install core packages
```bash
pip install numpy pandas scikit-learn matplotlib openpyxl imbalanced-learn
```


Download PreRisk-CoV2 scripts
```bash
git clone https://github.com/NTOUBiomedicalAILAB/PreRisk-CoV2.git
cd PreRisk-CoV2/
```

---

## Quick test
## Clinical datasets should follow ethical and privacy regulations
## Quick test
> Clinical datasets should follow ethical and privacy regulations

### Internal Validation (LOOCV)
```bash
python prerisk_cov2.py --mode internal --input Discovery.csv --n-iterations 100 --use-smote --plot-curves --output-dir ./results
```

### External Validation
```bash
python prerisk_cov2.py --mode external --train-input Discovery.csv --test-input Validation.csv --n-iterations 100 --use-smote --plot-curves --output-dir ./results
```


#### 2. Subsequent Usage

If the runs without errors, you only need to activate your environment before using PreRisk-CoV2:

```bash
conda activate PreRisk_CoV2
cd PreRisk-CoV2/
```

---

## 📊 Input Data Format
The input consists of protein expression data (CSV format),
To ensure compatibility with the prediction pipeline, please format your input CSV as follows:

### CSV File Structure

- **Column 0**: `sample ID` - Unique identifier for each patient/sample.
- **Column 1**: `PCR result` - Ground truth labels (can be `Detected`/`Not` or `1`/`0`).
  - *Note: If using `--no-labels` for pure prediction, this column can contain placeholders.*
- **Column 2 ~ N**: Protein expression levels (e.g., Olink NPX values).

### 🧬 The 7-Protein Panel (Default)
By default, the system automatically extracts the following 7 biomarkers using case-insensitive name matching:
> **MCP-3, LIF-R, TRANCE, FGF-23, NT-3, CXCL1, CXCL6**

<br>

---





## 📊 Data Availability

### Public Datasets

De-identified individual participant data supporting the findings of this study are available in the Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo) under accession numbers **GSE198449 (CHARM cohort)** and **GSE178967 (CEIM cohort)**.

1. NPX data for the **CHARM cohort** are accessible via the supplementary material of Soares-Schanoski et al. https://pmc.ncbi.nlm.nih.gov/articles/PMC9037090
2. NPX data for the **CEIM cohort** are available at the following repository: [https://github.com/hzc363/COVID19_system_immunology/blob/master/OLINK\%20Proteomics/olink.csv](https://github.com/hzc363/COVID19_system_immunology/blob/master/OLINK%20Proteomics/olink.csv)






---



