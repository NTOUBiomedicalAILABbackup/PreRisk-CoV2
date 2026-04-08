![PreRisk-CoV2](logo_v2.png)

## Overview

PreRisk-CoV2 is a machine learning framework for pre-exposure risk assessment of SARS-CoV-2 susceptibility using Serum protein biomarkers. The main function is to predict infection risk **before exposure** based on a 7-protein panel identified through K-Nearest Neighbors (KNN) combined with Genetic Algorithm (GA) feature selection. The input consists of protein expression data (CSV format), and the output provides risk prediction results with comprehensive performance metrics.

📄 **Paper**: Pre-exposure Serum Proteins Predict 3-Month SARS-CoV-2 Susceptibility: A 7-Protein Candidate Signature

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
pip install numpy pandas scikit-learn matplotlib openpyxl imbalanced-learn tenseal
```


Download PreRisk-CoV2 scripts
```bash
git clone https://github.com/NTOUBiomedicalAILABbackup/PreRisk-CoV2.git
cd PreRisk-CoV2/
```

---

## Quick test
## Clinical datasets should follow ethical and privacy regulations
Step 1: Encrypts your CSV to a secure .bin format
```bash
python encrypt.py --input Validation.csv --output Validation_encrypted.bin
```


Step 2: Set up encrypted-style KNN computation 

```bash
python prerisk_cov2.py --mode run --input Validation_encrypted.bin
```

Step 3: Decrypt the encrypted results upon request to generate the final report

- Step 3A: Upload your encrypted_result.bin to Google Drive
- Step 3B: Set sharing permissions  Right-click → Share → "Anyone with the link" → Viewer
- Step 3C: Copy the shareable link
- Step 3D: Send this exact email:
 ```txt
To: ntoubiomedicalai2026@gmail.com
Subject: [PreRisk-CoV2] Analysis Request
Body: [Paste your Google Drive link here]
```



Report Contains  

- Sample ID  
- Infection Risk Prediction (High/Low)  
- Risk Score (0.0-1.0)  


#### 2. Subsequent Usage

If the runs without errors, you only need to activate your environment before using PreRisk-CoV2:

```bash
conda activate PreRisk_CoV2
cd PreRisk-CoV2/
```

---


## Expected Response Time

**⏱️ Server will reply within 2 minutes** with:
- ✅ `Client_Report_Decrypted.csv` (analysis results)
- ✅ Privacy confirmation (your data deleted)



---


## 📊 Input Data Format
The input consists of protein expression data (CSV format), and the output provides infection risk prediction.
To ensure compatibility with the FHE encryption and prediction pipeline, please format your input CSV as follows:

### CSV File Structure

- **Column 0**: `sample ID` - Unique identifier for each patient/sample.
- **Column 1**: `PCR result` - Ground truth labels (can be `Detected`/`Not` or `1`/`0`).
  - *Note: If using `--no-labels` for pure prediction, this column can contain placeholders.*
- **Column 2 ~ N**: Protein expression levels (e.g., Olink NPX values).

### 🧬 The 7-Protein Panel (Default)
By default, the system automatically extracts and encrypts the following 7 biomarkers using  name matching:
> **MCP-3, LIF-R, TRANCE, FGF-23, NT-3, CXCL1, CXCL6**

<br>

---





## 📊 Data Availability

### Public Datasets

De-identified individual participant data supporting the findings of this study are available in the Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo) under accession numbers **GSE198449 (CHARM cohort)** and **GSE178967 (CEIM cohort)**.

1. NPX data for the **CHARM cohort** are accessible via the supplementary material of Soares-Schanoski et al. https://pmc.ncbi.nlm.nih.gov/articles/PMC9037090
2. NPX data for the **CEIM cohort** are available at the following repository: [https://github.com/hzc363/COVID19_system_immunology/blob/master/OLINK\%20Proteomics/olink.csv](https://github.com/hzc363/COVID19_system_immunology/blob/master/OLINK%20Proteomics/olink.csv)






---



