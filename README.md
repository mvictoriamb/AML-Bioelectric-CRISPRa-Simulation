# Leukemia Biomarker Discovery & Bioelectric Analysis
## 📝 Research Overview
This repository contains the computational pipeline and bioinformatics tools developed for the research paper: "Deep Learning and Bioelectric Analysis for Biomarker Discovery in Leukemia".
The project implements a Deep Learning framework for clinical subtype classification and feature selection. It specifically explores the synergy between genomic "dark matter" (lncRNAs) and cancer bioelectricity by identifying regulatory axes between non-coding transcripts and ion channels.

## 📂 Project Structure
├── data/
│   ├── Leukemia_GSE9476.csv       # Main dataset (Gene expression matrix)
│   └── GPL570-55999.txt           # Annotation file (Must be downloaded manually)
├── models/
│   └── modelo_leucemia.keras      # Trained Neural Network model
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
├── requirements.txt               # Required Python libraries
└── physicell/
    ├── output_control/
    ├── output_crispra/
    ├── output_rosi/
    └── plot_resultados.py

##📂 Data Requirements

Due to GitHub's file size restrictions, the platform annotation file is not included in this repository.

Download the file GPL570-55999.txt from NCBI GEO (GPL570).

Rename the file (if necessary) to GPL570-55999.txt.

Place it inside the /data folder in the project root.

This file is critical for mapping technical probe IDs to biological gene symbols.

*Note: The PhysiCell output folders contain the final state data and spatial SVG representations used to generate the figures in the manuscript. The full C++ custom modules are being refactored for a future standalone software release.*
