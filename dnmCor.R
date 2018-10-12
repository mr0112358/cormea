library(parallel)


dnmCor = function(simdata,nsamp = 199){
  
  
  n = length(simdata)
  g = length(simdata[[1]][[3]])-1
  
  # Extract Data
  X = t(vapply(simdata,function(x) x[[1]],numeric(g+1)))
  Y = t(vapply(simdata,function(x) x[[2]],numeric(g+1)))
  
  # Coefficients in Simpson's formula
  SimpVec = 1/(g*3)*c(1,4,rep(c(2,4),times = (g-2)/2),1)
  
  # Standardize X's and Y's
  Xcup = X - matrix(rep(colMeans(X),each = n),ncol = g+1)
  Mx = Xcup %*% SimpVec
  diffXcup = Xcup - matrix(rep(Mx,times=g+1),ncol = g+1)
  Xstar = matrix(rep(1/sqrt(diffXcup^2 %*% SimpVec),times=g+1),ncol = g+1) * diffXcup
  
  Ycup = Y - matrix(rep(colMeans(Y),each = n),ncol = g+1)
  My = Ycup %*% SimpVec
  diffYcup = Ycup - matrix(rep(My,times=g+1),ncol = g+1)
  Ystar = matrix(rep(1/sqrt(diffYcup^2 %*% SimpVec),times=g+1),ncol = g+1) * diffYcup
  
  # Here we use |rho| rather than rho because we just want to test dependency
  dnmcor = mean((Xstar * Ystar) %*% SimpVec)
  
  # Permutation test
  gmap = split(replicate(nsamp,sample.int(n)),rep(1:nsamp, each = n))
  alldnmcor = mclapply(gmap, function(perm){
    mean((Xstar * Ystar[perm,]) %*% SimpVec)
  },mc.cores = detectCores())
  pvalue = sum(abs(unlist(alldnmcor, use.names = FALSE))>=abs(dnmcor))/nsamp
  list(dnmcorr = dnmcor, pvalue = pvalue, unlist(alldnmcor, use.names = FALSE))
}

