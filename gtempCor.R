library(parallel)

gtempCor = function(simdata,nsamp = 199){


n = length(simdata)
g = length(simdata[[1]][[3]])-1

# Extract Data
X = t(vapply(simdata,function(x) x[[1]],numeric(g+1)))
Y = t(vapply(simdata,function(x) x[[2]],numeric(g+1)))

# Coefficients in Simpson's formula
SimpVec = 1/(g*3)*c(1,4,rep(c(2,4),times = (g-2)/2),1)

# Here we use |rhot| rather than rhot because we just want to test dependency
rhot = vapply(1:(g+1), function(j) cor(X[,j],Y[,j]), numeric(1))
rho = SimpVec %*% rhot

# Permutation test
gmap = split(replicate(nsamp,sample.int(n)),rep(1:nsamp, each = n))
allgtempcor = mclapply(gmap, function(perm){
  SimpVec %*% vapply(1:(g+1), function(j) cor(X[,j],Y[perm,j]), numeric(1))
},mc.cores = detectCores())
pvalue = sum(abs(unlist(allgtempcor, use.names = FALSE))>=abs(as.numeric(rho)))/nsamp
list(gtemprho = rho, pvalue = pvalue, unlist(allgtempcor, use.names = FALSE))
}