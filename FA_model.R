# factor_analysis function
#
# This function will compute the final factor analysis solution for a dataset combining all steps of the 
# recommended protocol into one function.
# 
# 1. choose correlation or covariance matrix
#
# 2. decide method of estimation (i.e. principal component, max likelihood, etc.) based on multivariate 
#    normality if instructed to do so
#
# 3. decide the number of factors to extract (VSS, PA, MAP, or 'systematic tuning')
#
# 4. decide type of rotation (i.e. oblique vs. varimax) based on outcome of correlation analysis
#
# 5. calculate final solution (output loadings to excel)
#
# 6. calculate factor scores (output scores to excel)
#
# 7. create output structure to save all necessary info for interpreation and reporting
#
#######################################################################################################################

# PARAMETERS (factor_analysis)
#
# X <- data matrix
# mat <- type of matrix ('cor' or 'cov')
# fm <- method of estimation for FA (will be used for final solution if check_normality set to FALSE)
#       options: 'pc', 'mle', 'pa', 'minres'
# check_normality <- use test of multivariate normality to choose method of estimation
# rotOrtho <- orthogonal rotation
# rotOblique <- oblique rotation
# extractMethod <- method to determine number of factors to extract ('sys','vss','map','parallel')
# iterations <- for systematic tuning bootstrapping procedure, unused if different extractMethod is chosen
# scores_method <- method of estimation for factor scores (see factor.scores in "psych" package 
# for full documentation and options))

factor_analysis <- function(X, mat='cor', fm = 'pc', check_normality = TRUE, rotOrtho='varimax',
                            rotOblique='promax',extractMethod = 'sys', iterations = 1000, 
                            scores_method='Thurstone'){
  
  # define dims
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # CHOOSE matrix for FA
  if (mat == 'cor'){
    faMat <- cor(X,use = "pairwise.complete.obs")
  } else if (mat == 'cov'){
    faMat <- cov(X,use = "pairwise.complete.obs")
  }
  
  # DETERMINE normality
  if (check_normality == TRUE){
    chiTest <- multivariate_normality(X)
    
    # DECIDE factor analysis model
    if (chiTest$p.value < 0.05){
      # non-parametric method: must use principal axis
      fm = 'pc'
      
    } else if (chiTest$p.value >= 0.05) {
      # parametric method: can change to whichever parametric method you'd like!
      fm = 'mle'
    }
  }
  
  # DETERMINE number of factors to extract
  if (extractMethod == 'sys'){
    
    sysNF = sysSolution(n,faMat,fm=fm,runLim=3,iterations=iterations)
    nF = sysNF[[1,1]]
    
  } else if (extractMethod == 'map'){
    nFAll <- vss(faMat,n=nF,rotate=rotate,fm=fm,n.obs=n, title="Number of Factors")
    nF <- which(nFAll$map == min(nFAll$map))
  } else if (extractMethod == 'vss'){
    nFAll <- vss(faMat,n=nF,rotate=rotate,fm=fm,n.obs=n, title="Number of Factors")
    nF <- which(nFAll$vss.stats$cfit.2 == max(nFAll$vss.stats$cfit.2))
  } else if (extractMethod == 'parallel'){
    nFpa2 = paran(X)
    nF = nFpa2$Retained
  }
  
  # CHOOSE rotation
  rotate_output = choose_rotation(faMat,n,fm,nF)
  rotate = rotate_output[[2,1]]
  
  #CALCULATE final model
  if (fm == 'pc'){
    
    if (rotate == rotOrtho){
      faFinal <- principal(faMat,nfactors=nF,residuals=TRUE,rotate=rotate,scores=FALSE,n.obs=n,eps=1e-14)
    } else if (rotate == rotOblique){
      faFinal <- principal(faMat,nfactors=nF,residuals=TRUE,rotate=rotate,n.obs=n,scores=FALSE)
    }
    
  } else {
    faFinal = fa(faMat,nF,n.obs=n,fm=fm,rotate=rotate,residuals=TRUE)
    
  }
  
  #calculate factor scores
  L = matrix(faFinal$loadings,nrow=p,ncol=nF)
  scores = factor.scores(X,L,method=scores_method)
  
  #output loadings and scores
  write.xlsx(L, "loadings.xlsx")
  write.xlsx(scores, "scores.xlsx")
  
  # save output
  output = matrix(list(),nrow=7,ncol=1)
  output[[1]] = fm
  output[[2]] = rotate
  output[[3]] = faFinal
  output[[4]] = chiTest
  output[[5]] = rotate_output[[1,1]]
  output[[6]] = scores$scores
  
  if (extractMethod == 'sys'){
    output[[7]] = sysNF
  } else {
    output[[7]] = extractMethod
  }
  
  output = data.frame(output,row.names = c('method','rotation','FA','normality','rotation corr','scores','extraction method'))
  
  return(output)
  
}

# OUTPUT structure describing the final and statistically validated solution
#
# output[[1,1]] = method of estimation
# output[[2,1]] = rotation
# output[[3,1]] = final model results
# output[[4,1]] = multivariate normality test
# output[[5,1]] = correlation results between oblique and orthogonal solutions
# output[[6,1]] = factor scores
# output[[7,1]] = method for choosing number of factors to extract (include systematic tuning solution if 
# this method is chosen)
#
##############################################################################################################

