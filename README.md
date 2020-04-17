# factor_analysis
## Factor analysis protocol

This repo includes code for a full recommended protocol to guide factor analyses. 

The faModel.R script will execute the following steps:

1. decide method of estimation (i.e. principal component, max likelihood, etc.)

2. decide the number of factors to extract (i.e. VSS, PA, MAP, boostrapping)

3. decide type of rotation (i.e. oblique vs. varimax)

4. compute a final factor model and visualize solution

This is accomplished using the following custom function libraries:

1. FA_utils.R - miscellaneous functions that are used throughout the factor analysis model
   (loadPacks, multivariate_normality, impute_data, choose_rotation, visualize_solution, factor_corrplot)
    
2. sysCriterion_utils.R - functions that quantify the systematic tuning criterion and the probability a
   given solution could have emerged by chance
   (sysFactors, sysSolution)
    
3. Fa_model.R - a function that pulls together all previously mentioned functions to implement the full 
   factor analysis protocol in one line of code
   (factor_analysis)
   
Datasets included:

1. hue-scaling
2. motion-scaling
