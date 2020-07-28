# KJE: This script is a protocol that should guide factor analyses generally, but specifically
# for scaling experiments. The script will execute the following steps:

# 1. decide method of estimation (i.e. principal component, max likelihood, etc.)

# 2. decide the number of factors to extract (i.e. VSS, PA, MAP, boostrapping)

# 3. decide type of rotation (i.e. oblique vs. varimax)

# Conclude with a final factor model and visualize solution

#######################################################################################################################
if ((!require(rstudioapi)) == TRUE){
  install.packages("rstudioapi")
} else {
  library(rstudioapi)
}

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

# LOAD function libraries
source('sysCriterion_utils.R')
source('FA_utils.R')
source('FA_model.R')

# CHECK for and install necessary packages
loadPacks()

######################################################################################################################

# LOAD and format data as X
filename <- 'finalData_motion.xlsx'
scalingData <- read_excel(filename,col_names = TRUE)
X <- matrix(as.numeric(unlist(scalingData)),nrow=dim(scalingData)[1])
faMat <- cor(X, use = "pairwise.complete.obs")

######################################################################################################################

# CHECK for correlation structure
corrplot(faMat,tl.col = 'black')

# BARTLETT'S test of sphericity is a quick check which tests the extent to which the correlation matrix resembles an 
# identity matrix. If it does, factor analysis is not the appropriate model  (pvalue < 0.05 indicates correlation 
# structure exists and that FA should be used)
sphericity <- cortest.bartlett(faMat,n=dim(X)[1])
sphericity$p.value

# KAISER MEYER OLKIN measure of samplng adequacy (rule of thumb, overall MSA > 0.5)
kmo_test <- KMO(faMat)
kmo_test$MSA

# if p.value is == NaN, check whether the det(R) is near 0. If so, the sphericity test can't be computed and the 
# determinant being near 0 suggests linearly dependent variables and thus supports the use of factor analysis.
det(faMat)

################################################################################################################

# CHECK for multivariate normality
multivariate_normality(X)

################################################################################################################

# CHECK scree plot
eig <- eigen(faMat)
plot(1:dim(faMat)[1],type='b',eig$values,main='scree plot',xlab='component',ylab='eigenvalue')

# RUN factor analysis
solution = factor_analysis(X, check_normality = FALSE, iterations=1000)

# VIEW model output
faFinal = solution[[3,1]]

################################################################################################################

#VISUALIZE solution
visualize_solution(scalingData,faFinal)

#VISUALIZE correlation plot for a single factor
factor_corrplot(scalingData,faFinal,1)

################################################################################################################

# CHECK solution(s) of more or less factors
nF = dim(faFinal$loadings)[2]
solution_plus = factor_analysis(X, extractMethod = 'none', nF = nF + 1, check_normality = FALSE, iterations=1000)
faFinal_plus = solution_plus[[3,1]]
visualize_solution(scalingData,faFinal_plus)


