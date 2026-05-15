# Leukemia Biomarker Discovery & Bioelectric Analysis
## 📝 Research Overview
This repository contains the computational pipeline and bioinformatics tools developed for the research paper: "In Silico Validation of CRISPRa-Mediated LOC401317 Activation to Reverse KCNJ15-Induced Bioelectric Immunosuppression in Acute Myeloid Leukemia".

The project implements a Deep Learning framework for clinical subtype classification and feature selection. It specifically explores the synergy between genomic "dark matter" (lncRNAs) and cancer bioelectricity by identifying regulatory axes between non-coding transcripts and ion channels.

This study validates the LOC401317-KCNJ15 axis as a critical bioelectric vulnerability in AML. By demonstrating that restoring the cellular membrane potential reverses tumor-induced immunosuppression, we provide a robust computational roadmap for therapeutic intervention. These findings establish a rigorous mechanistic foundation to transition from in silico modeling to experimental in vitro assays, highlighting a promising non-canonical pathway for leukemia treatment.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20208524.svg)](https://doi.org/10.5281/zenodo.20208524)

## 📂 Project Structure
````
├── data/
│   ├── Leukemia_GSE9476.csv       # Main dataset (Gene expression matrix)
│   └── GPL570-55999.txt           # Annotation file (Must be downloaded manually)
├── models/
│   └── modelo_leucemia.keras      # Trained Neural Network model
└── physicell/
│   ├── output_control/
│   ├── output_crispra/
│   ├── output_rosi/
│   └── plot_resultados.py
├── scripts/
│   ├── Red_Neuronal_CuMiDa.ipynb      # Main development notebook
│   ├── analisis_subtipos.py       # Generates expression heatmaps per clinical subtype
│   └── validacion_genomica.py      # Statistical validation (Boxplots & Pearson)
├── plots/
│   ├── correlacion_kcnj15_loc401317.png
|   ├── grafica_resultados_AML.png
│   ├── heatmap_subtipos_top20.png
│   ├── reporte_interactivo.html    # Interactive Plotly report
│   ├── sobreexpresion_kcnj15.png
│   └── sobreexpresion_loc401317.png
├── results/
│   └── biomarcadores.csv          # Top genes identified by the AI 
├── CITATION.cff
├── LICENSE
├── README.md
└── requirements.txt               # Required Python libraries
````

## 📂 Data Requirements

Due to GitHub's file size restrictions, the platform annotation file is not included in this repository.

Download the file GPL570-55999.txt from NCBI GEO (GPL570).

Rename the file (if necessary) to GPL570-55999.txt.

Place it inside the /data folder in the project root.

This file is critical for mapping technical probe IDs to biological gene symbols.

*Note: The PhysiCell output folders contain the final state data and spatial SVG representations used to generate the figures in the manuscript. The full C++ custom modules are being refactored for a future standalone software release.*
