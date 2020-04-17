############################################################################################################
# SYSTEMATIC TUNING CRITERION UTILS
# Set of functions used to quantify the systematic tuning criterion. Specifically, these functions run 
# the bootstrapping analyses to determine that probability that a factor with a given number of consecutive 
# loadings over a set threshold could have emerged by chance. An algorithm for finding a stastically 
# significant solution for a given dataset is also included.
#
# 1. sysFactors <- runs the bootstrapping analysis to determine that probability that a factor with a given 
# number of consecutive loadings over a set threshold could have emerged by chance
#
# 2. sysSolution <- algorithm that finds a stastically significant solution for a given dataset
# 
#############################################################################################################

# PARAMETERS (sysFactors): 
#
# n <- # of samples 
# p <- # of variables 
# fm <- method of estimation for factor model 
# rotate <- type of rotation for factor model
# nF <- # of factors to extract
# ldgTh <- loading threshold
# runTh < # of consecutive loadings 
# iterations <- # of iteration for bootstrapping 
# 
# NOTE: (n,p,fm,rotate) should match the original dataset as closely as possible

sysFactors <- function(n,p,fm,rotate,nF,ldgTh,runTh,iterations=1000){

    allSys = as.numeric()
  
  for (i in 1:iterations){
    
    # correlation matrix of random data
    corMat <- cor(matrix(rnorm(n*p),n),use = "pairwise.complete.obs")
    
    # factor analysis
    if (fm == 'pc'){
      factorSolution <- principal(corMat,nfactors=nF,residuals=TRUE,rotate=rotate,n.obs=n,eps=1e-14)
      
    } else {
      factorSolution = fa(corMat,nF,n.obs=n,fm=fm,rotate=rotate,residuals=TRUE,scores="regression")
      
    }
    
    sig_loadings = matrix(as.numeric(abs(factorSolution$loadings) > ldgTh),ncol=nF)
    
    #take a count of the number of 'systematic factors'
    sysFactors = as.numeric()
    
    for (factor in 1:nF){
      
      #find runs of loadings greater than threshold
      runs = rle(sig_loadings[,factor])
      vals = as.numeric(which((runs$values == 1) == TRUE))
      lens = as.numeric(which((runs$lengths >= runTh) == TRUE))
      ans = sum(as.numeric(vals %in% lens))
      
      #check runs
      if (ans > 0){
        runIDX = lens[lens %in% vals]
        initIDX = sum(runs$lengths[1:runIDX-1])+1
        runLoadings = factorSolution$loadings[initIDX:(initIDX+runs$lengths[runIDX]-1),factor]
        
        #redo calc of runs based on sign of loadings
        sign_runs = rle(sign(runLoadings))
        sr = sum(sign_runs$lengths >= runTh)
        
        if (sr > 0){
          sysFactors[factor] = 1
        } else {
          sysFactors[factor] = 0
        }
      } else{
        sysFactors[factor] = 0
      }
    }
    
    allSys[i] = sum(sysFactors)
  }
  
  # calculate and return p-value for the given loading and consecutive loading thresholds  
  sysFreq = sum(allSys)/(iterations*nF)
  return(sysFreq)
}

######################################################################################################################

# PARAMETERS (sysSolution): 
#
# n <- # of samples 
# faMat <- data correlation/covariance matrix
# fm <- method of estimation for factor model 
# rotate <- type of rotation for factor model
# initNF <- initial number of factors (Kaiser rule unless otherwise specified)
# runLim <- limit for number of consecutive loadings to search for in solution
# iterations <- # of iteration for bootstrapping 
# 
# NOTE: (n,fm,rotate) should match the original solution

sysSolution <- function(n,faMat,fm,rotate="varimax",initF = 0,runLim=4,iterations=1000){
  
  #unless otherwise specified, nF will be initialized as the number of eigenvalues greater than 1 
  # (aka greater than average)
  if (initF == 0){
    eigs <- eigen(faMat)
    nLambda <- sum(eigs$values > 1)
    nF <- nLambda
  } else {
    nF = initF
  }
  
  
  runTh = 2
  pSys = 1
  
  while (pSys >= 0.05){
    
    #compute factor model with given parameters
    if (fm == 'pc'){
      
      if (rotate == "varimax"){
        faNF <- principal(faMat,nfactors=nF,residuals=TRUE,rotate=rotate,scores=FALSE,n.obs=n,eps=1e-14)
      } else {
        faNF <- principal(faMat,nfactors=nF,residuals=TRUE,rotate=rotate,n.obs=n,scores=FALSE)
      }
      
    } else if (fm == fmPara){
      faNF = fa(faMat,nF,n.obs=n,fm=fm,rotate=rotate,residuals=TRUE)
      
    }
    
    #find highest systematic loadings per factor
    sysLoadings = matrix(list(),nrow=nF,ncol=1)
    
    for (factor in 1:nF) {
      ldgs = faNF$loadings[,factor]
      maxn = sort(abs(ldgs),decreasing=TRUE)
      idx = 1
      
      while(sum(sysLoadings[[factor]]) == 0){
        max_loadings = (ldgs >= maxn[idx]) + (ldgs <= -maxn[idx])
        runs = rle(max_loadings)
        
        #check if high loadings are part of systematic sequence
        if (runs$lengths[which(runs$values == 1)] >= runTh){
          cumIDX = cumsum(runs$lengths)
          runIDX = which(runs$values == 1)
          sysLoadings[[factor]] = ldgs[(cumIDX[runIDX-1]+1):cumIDX[runIDX]]
        } else {
          idx = idx + 1
        }
        
      }
    }
    
    #evaluate the chance probability of the loading and run thresholds that characterize the current model
    #as long as runs for all factors exist, otherwise automatically search for different model
    if (sum(which(lengths(sysLoadings) < runTh)) == 0) {
      ldgTh = min(unlist(sysLoadings))
      pSys = sysFactors(dim(X)[1],dim(X)[2],fm,rotate,nF,ldgTh,runTh,iterations)
    } else {
      pSys = 1
    }
    
    #if this prob is too high, add a run until limit, then decrease number of factors
    if (pSys > 0.05){
      runTh = runTh + 1
      if (runTh > runLim){
        runTh = 2
        nF = nF - 1
      }
    }
    
  }
  
  output = matrix(list(),nrow=5,ncol=1) 
  output[[1]] = round(nF)
  output[[2]] = pSys
  output[[3]] = ldgTh
  output[[4]] = round(runTh)
  output[[5]] = sysLoadings
  output = data.frame(output,row.names=c('# of factors', 'pval', 'loading thresh','consecutive loadings','sys loadings'))
  
  return(output)
}

# OUTPUT structure describing the final and statistically significant solution
#
# output[[1,1]] = # of factors
# output[[2,1]] = p-value 
# output[[3,1]] = loading threshold
# output[[4,1]] = consecutive loading threshold
# output[[5,1]] = loadings for each factor that met the systematic tuning criterion
#
##############################################################################################################