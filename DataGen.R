## Loading packages
library(mvtnorm)

#Covariance Structure #Suppose m=8, k=1,...,m
m=8
Rho = diag(c(0,0,0,0,0,0,0,0))
betaX = 2
betaY = 3
Sigma11 = diag(4*seq(1,m)^(-betaX))
Sigma22 = diag(6*seq(1,m)^(-betaY))
Sigma12 = sqrt(Sigma11) %*% Rho %*% sqrt(Sigma22)
Sigma = rbind(
  cbind(Sigma11, Sigma12),
  cbind(t(Sigma12), Sigma22)
)


# Generating \xi's and \zeta's for n times, n=100 or 500 or ...
n=100
Coef = rmvnorm(n,sigma = Sigma)
xi = Coef[,1:m]
xi = lapply(1:n,function(x) t(xi[x,]))
zeta = Coef[,(m+1):(2*m)]
zeta = lapply(1:n,function(x) t(zeta[x,]))
rm(Coef)

# Basis functions
sqrt2 = sqrt(2)
X_basis = function(t){
  sqrt2*rbind(cos(2*pi*t),
              cos(4*pi*t),
              cos(6*pi*t),
              cos(8*pi*t),
              cos(10*pi*t),
              cos(12*pi*t),
              cos(14*pi*t),
              cos(16*pi*t))
}
Y_basis = function(t){
  sqrt2*rbind(
    sin(2*pi*(t+.2)),
    sin(4*pi*(t+.2)),
    sin(6*pi*(t+.2)),
    sin(8*pi*(t+.2)),
    sin(10*pi*(t+.2)),
    sin(12*pi*(t+.2)),
    sin(14*pi*(t+.2)),
    sin(16*pi*(t+.2))
  )
}

# Sampling Times
g = 40
Time = rep(list((0:g)/g),times=n)

# Simulation Data Generating function
simDataGen = function(basisX,basisY,coefsX,coefsY,Time){
  ##############################################################################
  # Input:
  ##    basisX and basisY are basis functions of X's and Y's respectively
  ##    coefsX and coefsY are coefficients of basis function of X's and Y's respectively
  ##    time is a list of sampling time for subject i=1,...,n
  
  # Output:
  ## A list of n sublist which contains X_ij, Y_ij and time t_ij,
  
  ##############################################################################
  mapply(
    FUN = list,
    lapply(mapply(FUN=list,lapply(Time,basisX),coefsX,SIMPLIFY = FALSE), 
           function(x) {as.vector(x[[2]] %*% x[[1]])}),
    lapply(mapply(FUN=list,lapply(Time,basisY),coefsY,SIMPLIFY = FALSE), 
           function(x) {as.vector(x[[2]] %*% x[[1]])}),
    Time,
    SIMPLIFY = FALSE
  )
}

simdata = simDataGen(X_basis,Y_basis,xi,zeta,Time)

rm(list = (ls()[ls()!="simdata"]))
