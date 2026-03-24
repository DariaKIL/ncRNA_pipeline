# ncRNA_pipeline

Computational pipeline for the analysis and visualization of small non-coding RNA (microRNA and tRNA-derived fragments) in diverse biological datasets.

This repository contains reusable analysis functions, visualization utilities, and workflow structure designed to process multiple ncRNA datasets generated in cardiovascular and transplantation-related studies.

The pipeline is being developed to support several ongoing research projects involving extracellular vesicles, heart transplantation cohorts, and pancreatic cyst-derived biological material.

⚠️ **Status: under active development**
The repository is intended as a general analysis framework. Individual research projects and publications may have their own dedicated repositories with finalized analysis scripts and curated datasets.

*Project overview*

This repository contains analysis workflows and exploratory data analysis for several datasets involving non-coding RNA expression:
- extracellular vesicles from cardiovascular samples
- heart transplant recipient cohorts
- pancreas cyst-derived biological material

[Data visualization transplantations](https://dariakil.github.io/ncRNA_pipeline/Transplantation/reports/Transplantation_analysis.nb.html)

**Repository structure (planned)**
```
ncRNA_pipeline/
├── README.md
├── .gitignore
├── .gitattributes
├── config.R
├── data/
│   └── raw/
│       ├── miR.Counts.csv
│       ├── tRF.Counts.csv
│       └── annotation.report.csv
├── src/
│   ├── data_processing.R
│   ├── deseq_analysis.R
│   ├── enrichment.R
│   ├── visualization.R
│   ├── generate_figures.R
│   ├── wgcna_analysis.R
│   └── calculations.R
├── Transplantation/
│   ├── config_t.R
│   ├── data/
│   │   ├── phenotable.tsv
│   │   └── Cytokines.csv
│   ├── Transplantation_analysis.Rmd
│   ├── Transplantation_WGCNA.Rmd
│   ├── reports/
│   │   ├── Transplantation_analysis.nb.html
│   │   └── Transplantation_WGCNA.nb.html
│   └── figures/
├── Cysts_pancreas/
│   ├── config_c.R
│   ├── data/phenotable.tsv
│   ├── Cysts_analysis.Rmd
│   ├── reports/
│   │   └── Cysts_analysis.nb.html
│   └── figures/
└── Vesicles_CABG/
    ├── config_v.R
    ├── data/phenotable.tsv
    ├── Vesicles_analysis.Rmd
    ├── reports/
    │   └── Vesicles_analysis.nb.html
    └── figures/
``` 
