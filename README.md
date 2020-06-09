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
- Data contains all raw and processed data files used to run the analyses. Relevant data files are read in at the top of each script.
- Env data contains all raw, unprocessed environmental data. These data are processed in BI_abiotic_script.
- BI_RO_manip processes focal species based on their range occupancies.

## Analyses
- brms_BI_mod_LL contains the R files and R output used for the Bayesian modelling component of this project. 
- Figures contains each figure used in the manuscript in pdf format, as well as a powerpoint file used to make some manual aesthetic edits for the manuscript.
- Scripts includes all scripts used to complete this project. 
  - Data cleaning: 
    - BBS_occ_script downloads and processes the raw BBS data for the project
    - BI abiotic script processes the raw environmental data for the project
    - BI GIS script calculates the amount of overlap between each focal species and its competitor, which is used to select relevant competitors as well as identify the main competitor. Additionally, we calculate range size in this script.
    - Biotic_Interactions_Snell_Data_cleaning cleans and processes the BBS data, matches up the focal species with relevant competitors, and combines the BBS data with corresponding environmnetal data based on latitude and longitude.
    
  - Data analysis: 
    - Biotic_Interactions_Analysis (all comps) runs all analyses in the project other than the Bayesian models using the summed abundance of any competitor with range overlap with the focal species as the "competitor" variable.
    - Biotic_Interactions_Analysis (main comp) runs all analyses in the project other than the Bayesian models using the summed abundance of a single, main competitor with the most range overlap with the focal species as the "competitor" variable.
    - Biotic_Interactions_Analysis (post-hoc comps) runs all analyses in the project other than the Bayesian models using the summed abundance of any competitor identified as relevant in post-hoc anlayses with the focal species as the "competitor" variable.
    - phylogeny code runs a phylogenetic query and generates a tree using the focal species in our analysis.
  
