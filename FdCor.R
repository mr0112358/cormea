library(energy)
library(parallel)

FdCor = function(simdata,nsamp = 199){
  ##############################################################################
  # Input:
  ## simdata     Raw Data, a list of sublist contains X's, Y's and corresponding mesure times
  ## nsamp       the size for permutation test in a Monte Carlo method, default value = 199
  
  # Output:
  ## p-value
  ## All of FdCors by permutation
  
  # Note: in the FdCor, we require the sampling times for each subject are equal
  ## In this version 0.1, we assume sampling on grid time points on [0,1], i.e.,
  ## 0/g, 1/g, 2/g, ... , g/g, and g is a even number, so that we can use the 
  ## Simpson's numerical integration formula.
  ##############################################################################
  n = length(simdata)
  g = length(simdata[[1]][[3]])-1
  
  # Extract X's and Y's from simdata
  extX = function(j) {
    ## Input: time index j from 0 to g
    ## Output: A matrix, each row is the Xi(j)
    vapply(simdata,function(x) x[[1]][j] ,numeric(1))
  }
  
  extY = function(j) {
    ## Input: time index j from 0 to g
    ## Output: A matrix, each row is the Yi(j)
    vapply(simdata,function(x) x[[2]][j] ,numeric(1))
  }
  
  # Create coefficients matrix of Simpson's formula
  SimpVec = 1/(g*3)*c(1,4,rep(c(2,4),times = (g-2)/2),1)
  SimpMat = SimpVec %*% t(SimpVec)
  
  # calculate dcovXY, dvarX, dvarY for any time combination
  dcovXY = matrix(mapply(function(x,i,j){dcov(extX(i),extY(j))}, 
                         SimpMat,1:(g+1), 1:(g+1)),nrow=(g+1))
  dvarX  = matrix(mapply(function(x,i,j){dcov(extX(i),extX(j))}, 
                         SimpMat,1:(g+1), 1:(g+1)),nrow=(g+1))
  dvarY  = matrix(mapply(function(x,i,j){dcov(extY(i),extY(j))}, 
                         SimpMat,1:(g+1), 1:(g+1)),nrow=(g+1))
  
  # calculate numerical integration by Simpson's formula and get fdvarX, fdvarY and fdcovXY
  fdvarX  = sum(SimpMat*dvarX)
  fdvarY  = sum(SimpMat*dvarY)
  denom = sqrt(fdvarX*fdvarY)
  
  fdcovXY = sum(SimpMat*dcovXY)
  fdcorXY = fdcovXY/denom
  
  # permutation test, gmap is nsamp of permutation's of 1:n
  gmap = split(replicate(nsamp,sample.int(n)),rep(1:nsamp, each = n))
  allFdcor = mclapply(gmap, function(perm){
    dcovXY_perm = matrix(mapply(function(x,i,j){dcov(extX(i),extY(j)[perm])}, 
                                SimpMat,1:(g+1), 1:(g+1)),nrow=(g+1))
    sum(SimpMat*dcovXY_perm)/denom
  },mc.cores = detectCores())
  pvalue = sum(unlist(allFdcor, use.names = FALSE)>=fdcorXY)/nsamp
#  list(fdcor = fdcorXY, pvalue = pvalue,unlist(allFdcor, use.names = FALSE))
}