# Interplay between chromosomal alterations and gene mutations in clonal hematopoiesis
Code and Data for [Gao et al. Nature Communications 2021](https://www.nature.com/articles/s41467-020-20565-7).

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
|   |   └─ data/: example input data and commandline script for FACETS-CH
```
# Analysis notebooks
An interactive jupyter notebook is provided to reproduce the analyses in the paper. First install [Jupyter](https://jupyter.org/) and then launch the notebook [ch_cnv_figures.ipynb](https://github.com/papaemmelab/Gao_NC_CH/blob/main/notebooks/ch_cnv_figures.ipynb).
## Correction (03/22/2021)
The confidence level of confidence intervals (as specified in the figure legends) in Figure 1f, Figure 2b-c,  Supplementary Figure 9b, and Supplementary Figure 13 are 90% instead of 95%. No conclusions in the main or supplementary text are affected. 


# FACETS-CH algorithm
Code to implement the FACETS-CH algorithm is provided here under the [facets-ch](https://github.com/papaemmelab/Gao_NC_CH/tree/main/facets-ch) subfolder.
