# Biotic-Interactions
Determining whether biotic or abiotic factors more strongly determine avian species temporal occupancy, published in Global Ecology and Biogeography, 2020 (https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.13064)

The Biotic-Interactions repository contains five folders and an RProj file where analyses were run:
- brms_BI_mod_LL
- Data
- Env data
- Figures
- Scripts
- Biotic-Interactions.RProj

Throughout, "all" refers to scripts and data where the sum of all competitor abundances co-occurring with the focal species are used in the analysis, while "main" used a single competitor's abundance, chosen based on the amount of range overlap between the focal and competitor species.

## Data
The data folder contains all raw and processed data files used to run the analyses. Relevant data files are read in at the top of each script.
The env data folder contains all raw, unprocessed environmental data. These data are processed in BI_abiotic_script.

## Analyses
The brms_BI_mod_LL folder contains the R files and R output used for the Bayesian modelling component of this project. 
