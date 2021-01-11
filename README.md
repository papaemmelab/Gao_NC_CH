# Interplay between chromosomal alterations and gene mutations in clonal hematopoiesis
Code and Data for Gao et al. Nature Communications 2021.

![image](/onco.png)

# Directory structure
```
└─ notebooks/:
|   └─ ch_cnv_figures.html: Analysis code and scripts to reproduce figures
|   └─ plot_forest.R: helper script to create forest plots
|   └─ table_one.R: helper script to tally demographics
└─ data/:
|   └─ clinical_data.tsv: Patient blood count, demographics, and treatment data
└─ facets-ch/:
|   └─ R/: Source code for the FACETS-CH algorithm
|   └─ tests/:
|   |   └─ data/: example input data and commandline script for facets-ch
```
# Analysis notebooks
An interative jupyter notebook is provided to reproduce the analyses in the paper. First install [Jupyter](https://jupyter.org/) and then launch the notebook `ch_cnv_figures.ipynb`.
