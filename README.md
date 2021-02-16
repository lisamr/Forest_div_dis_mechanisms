# Forest_div_dis_mechanisms
Model code and data used for the manuscript: "Direction and underlying mechanisms of the diversity-disease relationship are distinct across hierarchical scales"

## Data is located in the "data" folder
- plot_level_data.csv = plot-level attributes  
- tree_level_data_HS.csv = tree-level attributes of highly/commonly symptomatic species only  
- tree_level_data_all.csv = tree-level attributes of all species  
- PA_matrix.csv = Presence/absence matrix of species  
- sporangia_estimates.csv = mean, sd of sporulation estimates; derived from Rosenthal et al. (in press, Plant Disease)  
- dmat.RDS = matrix of pairwise distances between plots, used for spatial disease risk models (see Appendix S1)  

Columns use accronyms or names as described in the manuscript.  

## Scripts 
Scripts used to run the models are located in the "scripts_models" folder and for figures/tables, in the "scripts_fig_tables". Stan models are located in the "Stan_models" folder.  

Analysis covers models used to assess:  

1. the relationship between various forms of density and diversity,   
2. community competence and diversity,  
3. factors that control disease risk at the individual and plot level. Included are he disease risk models used in the main text, as well as Gaussian process models, which were used to test if including a spatially weighted term changed the posteriors (it didn't).  
Files are named accordingly.  

I would run the model scripts first, save the model fits in another folder outside of this repository, and then read in the models fits again when creating the figures. Make sure to alter your directory when/if you want to save the models.  



----
Please contact Lrosenthal@ucdavis.edu for more information.  
