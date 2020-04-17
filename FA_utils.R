############################################################################################################
# FACTOR ANALYSIS UTILS 
# Library of function used across the factor analysis protocol
#
# 1. loadPacks <- checks whether to install and then loads all packages associated with the factor
#    analysis protocol
#
# 2. multivarate_normality <- tests whether the statistical distances of the dataset are chi-square 
#    distributed, as they should be for multivarate normal data
#
# 3. impute_data <- impute missing values by either non-parametric or parametric methods (depending 
#    on multivariate normality assumptions)
# 
# 4. choose_rotation <- choose a rotation for the final model based on the extent to which the Oblique 
#    and Orthogonal solutions are correlated. 
#
# 5. visualize_solution <- plot correlation matrices (original, model, residual) and factors
#
# 6. factor_corrplot <- plot individual factor and variance structure it accounts for
#
#############################################################################################################

# CHECK for packages to install
loadPacks <- function(){
  necessary_packs <- c("devtools","ggcorrplot","grid","gridExtra","paran","missForest","mice","foreach",
                       "doParallel","readxl","corrplot","ggplot2","psych","MASS","boot","lattice","nFactors",
                       "paran","matrixStats","plotly","ggplotify","GPArotation","openxlsx","easypackages",
                       "GPArotation")
  need2install <- necessary_packs %in% installed.packages() == FALSE
  
  if (sum(need2install) > 0) {
    install.packages(necessary_packs[need2install])
  }
  
  # load all necessary packages
  library(easypackages)
  
  libraries(necessary_packs)
  
  return()
}

# CHECK for multivariate normality (statiscal distances are chi-squared distributed)
#
# X <- data matrix
multivariate_normality <- function(X){
  
  d2 <- rep(0,dim(X)[1])
  S_1 <- solve(cov(X, use = "pairwise.complete.obs"))
  mu <- colMeans(X, na.rm = TRUE)
  
  for (i in 1:dim(X)[1]){
    d2[i] <- t(X[i,] - mu) %*% S_1 %*% (X[i,]-mu)
  }
  
  chiTest <- ks.test(d2,"pchisq",dim(X)[2])
  return(chiTest)
}


# IMPUTE missing values (check for normality and use parametric or non-parametric method accordingly)
#
# X <- data matrix
# xl_filename <- excel filename for imputed data
impute_data <- function(X,xl_filename="imputed_data.xlsx"){
  
  #check for normality
  chiTest = multivariate_normality(X)
  
  if (chiTest$p.value < 0.05){
    
    #non-parametric method
    imputed_data <- missForest(X)
    X <- imputed_data$ximp
    
  } else if (chiTest$p.value >= 0.05) {
    
    #parametric method
    imputed_data <- mice(X, m = 5, method = "pmm")
    X <- complete(imputed_data)
    X <- matrix(as.numeric(unlist(X)),nrow=n)
  }
  
  # SAVE data to excel for record of values that will be used for final analysis
  finalX = rbind(vars,X)
  write.xlsx(finalX, xl_filename)
  return(X)
}

# CHOOSE a rotation for the final model based on the extent to which the Oblique and Orthogonal solutions 
# are correlated. If they are correlated, this indicates an underlying orthogonal structure and therefore
# an orthogonal rotation is chosen.
#
# faMat <- data correlation/covariance matrix
# n <- # of samples
# fm <- method of estimation for factor analysis
# nF <- # of factors 
# rotOblique <- oblique rotation method
# rotOrtho <- orthogonal rotation method
#
# NOTE: (n, fm, nF) should match those of the solution for the original data
choose_rotation <- function(faMat,n,fm,nF,rotOblique='promax',rotOrtho='varimax'){
  
  p <- dim(faMat)[1]
  
  #compute factor analysis
  if (fm == 'pc'){
    faOblique <- principal(faMat,nfactors=nF,residuals=TRUE,rotate=rotOblique,n.obs=n,scores=FALSE)
    faOrtho <- principal(faMat,nfactors=nF,residuals=TRUE,rotate=rotOrtho,n.obs=n,eps=1e-14,scores=FALSE)
  } else {
    faOblique <- fa(faMat,nfactors=nF,n.obs=n,fm=fm,rotate=rotOblique,residuals=TRUE,scores="regression")
    faOrtho <- fa(faMat,nfactors=nF,n.obs=n,fm=fm,rotate=rotOrtho,residuals=TRUE,scores="regression") 
  }
  
  LOblique <- matrix(faOblique$loadings,ncol = nF*p)
  LOrtho <- matrix(faOrtho$loadings,ncol = nF*p)
  kstestOblique <- ks.test(LOblique,"pnorm",p)
  kstestOrtho <- ks.test(LOrtho,"pnorm",p)
  
  #correlate solutions 
  if (kstestOblique$p.value < 0.05 & kstestOrtho$p.value < 0.05){
    #non-parametric method
    corrRotate <- cor.test(LOblique,LOrtho,method="spearman")
  } else {
    #parametric method
    corrRotate <- cor.test(LOblique,LOrtho,method="pearson")
  }
  
  #if they are substantially correlated, then decide on orthogonal (e.g. varimax) rotation
  if (corrRotate$p.value < 0.05){
    rotate = rotOrtho
  }  else if (corrRotate$p.value >= 0.05){
    rotate = rotOblique #other oblique rotations should be tried, e.g. oblimin
  }
  
  output = matrix(list(),nrow=2,ncol=1)
  output[[1]] = corrRotate
  output[[2]] = rotate
  output = data.frame(output)
  
  return(output)
}


#VISUALIZE solution (correlation matrices: original, model, residual; factor plots)
#
# df <- data as data.frame
# faFinal <- final factor analysis solution from the factor_analysis function
visualize_solution <- function(df,faFinal){
  
  vars = as.numeric(colnames(df))
  
  #compute loading and uniquenesses matrices to calculate model estimated covariance structre
  L = matrix(faFinal$loadings,nrow=length(vars),ncol=dim(faFinal$loadings)[2])
  U = diag(faFinal$uniquenesses)
  nF = dim(L)[2]
  
  faMat <- cor(df)
  estMat <- L%*%t(L) + U
  colnames(estMat) <- vars
  rownames(estMat) <- vars
  resMat <- faMat - estMat
  
  L = data.frame(L)
  L$vars = vars
  
  head(faMat[, 1:length(vars)])
  p1 <- ggcorrplot(faMat,title='Data Correlation Matrix')
  
  head(estMat[,1:length(vars)])
  p2 <- ggcorrplot(estMat, title='Model Correlation Matrix')
  
  head(resMat[, 1:length(vars)])
  p3 <- ggcorrplot(resMat,title='Residual Correlation Matrix')
  
  head(estMat[,1:length(vars)])
  p2 <- ggcorrplot(estMat, title='Model Correlation Matrix')
  
  gs <- lapply(1:nF, function(i) 
    as.grob(ggplot(data=L, aes(x=vars, y=L[,i], group=1)) +
              geom_line()+
              geom_point()+
              labs(y = 'loading')+
              ylim(-1,1)))
  
  #windows()
  grid.arrange(grobs = list(p1,p2,p3),ncol =2)
  #windows()
  grid.arrange(grobs=gs, ncol=ceiling(nF/3))
}

#VISUALIZE factor plot and correlation structure accounted for by an individual factor
#
# df <- data as data.frame
# faFinal <- final factor analysis solution from the factor_analysis function
# i <- the factor to be visualized 
factor_corrplot <- function(df,faFinal,i){
  
  vars = as.numeric(colnames(df))
  L = matrix(faFinal$loadings,nrow=length(vars),ncol=dim(faFinal$loadings)[2])
  U = diag(faFinal$uniquenesses)
  
  Uind <- matrix(0, nrow=dim(U)[1], ncol=dim(U)[2])
  Uind[i,i] = U[i,i]
  estMat = L[,i]%*%t(L[,i]) + Uind
  colnames(estMat) <- vars
  rownames(estMat) <- vars
  
  L = data.frame(L)
  L$vars = vars
  
  p1 = ggplot(data=L, aes(x=vars, y=L[,i], group=1)) +
    geom_line()+
    geom_point()+
    labs(y = 'loading')+
    ylim(-1,1)
    
  head(estMat[,1:length(vars)])
  p2 = ggcorrplot(estMat, title=paste('Factor ',i,' Correlation Matrix'))
  
  grid.arrange(grobs = list(p1,p2),ncol =2,heights=c(1,1),widths=c(1,1))
}

