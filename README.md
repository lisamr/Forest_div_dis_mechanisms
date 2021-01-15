# Forest_div_dis_mechanisms
Model code and data used for the manuscript: "The importance of diversity and its inferred mechanisms on forest disease risk depend on how disease is measured"

## Data is located in the "data" folder
- plot_level_data.csv = plot-level attributes  
- tree_level_data_HS.csv = tree-level attributes of highly susceptible species only  
- tree_level_data_all.csv = tree-level attributes of all species   
- dmat.RDS = pairwise distance (km) matrix between plots (additional information for running a model that accounts for spatial autocorrelation. It wasn't necessary for this analysis)  

Columns use accronyms or names as described in the manuscript. 

## Scripts 
Analysis covers models used to assess   
1. the relationship between various forms of density and diversity and   
2. factors that control disease risk at the individual and plot level. Files are named accordingly.  

Scripts used to run the models are located in the "scripts_models" folder and for figures/tables, in the "scripts_fig_tables". Stan models are located in the "Stan_models" folder.  

----
Please contact Lrosenthal@ucdavis.edu for more information.  
