library(energy)
library(parallel)

FPCAdCor = function(simdata,nsamp = 199){
  
  n = length(simdata)
  g = length(simdata[[1]][[3]])-1
  
  # Extract X and Y data matrices
  X = t(vapply(simdata,function(x) x[[1]],numeric(g+1)))
  Y = t(vapply(simdata,function(x) x[[2]],numeric(g+1)))
  
  # Do PCA on X and Y respectively, then project centered X and Y on first n PCs 
  # such that the cumulant percentage is over 99
  prinX = prcomp(X)
  Mx = which.min(cumsum(prinX$sdev)/sum(prinX$sdev)<0.99)
  Xhat = prinX$x[,1:Mx]
  
  prinY = prcomp(Y)
  My = which.min(cumsum(prinY$sdev)/sum(prinY$sdev)<0.99)
  Yhat = prinY$x[,1:My]
  
# apply dcov in `energy` on Xhat and Yhat
  dvarX = dcov(Xhat, Xhat)
  dvarY = dcov(Yhat, Yhat)
  denom = sqrt(dvarX * dvarY)
  
  dcorXY = dcov(Xhat,Yhat)/denom
  
  # Permutation test
  gmap = split(replicate(nsamp,sample.int(n)),rep(1:nsamp, each = n))
  alldcor = mclapply(gmap, function(perm){
    dcov(Xhat,Yhat[perm,])/denom
  },mc.cores = detectCores())
  pvalue = sum(unlist(alldcor, use.names = FALSE)>=dcorXY)/nsamp
  #list(dcor = dcorXY, pvalue = pvalue, unlist(alldcor, use.names = FALSE))
}