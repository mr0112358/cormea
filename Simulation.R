# Load source files
source("./datagen.R")
source("./FdCor.R")
source("./FPCAdCor.R")
source("./dnmCor.R")
source("./gtempCor.R")

# number of observations
n = 100

# Sampling Times
g = 40
Time = rep(list((0:g)/g),times=n)

# Coeff parameters
Rho = diag(c(0,0,0,0,0,0,0,0))
betaX = 1.05
betaY = 1.2

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


# Coef = coefGen(Rho,betaX,betaY,n)
# simdata = simDataGen(X_basis,Y_basis,Coef,Time)

allsimdata = lapply(1:100,
       function(x){
         set.seed(x)
         Coef = coefGen(Rho,betaX,betaY,n)
         simdata = simDataGen(X_basis,Y_basis,Coef,Time)
       })

gtemp.out = vapply(allsimdata, gtempCor,numeric(1))
dnm.out = vapply(allsimdata,dnmCor,numeric(1))
FPCA.out = vapply(allsimdata,FPCAdCor,numeric(1))
FdCor.out = vapply(allsimdata,FdCor,numeric(1))
print(sum(gtemp.out<0.05)/100)
print(sum(dnm.out<0.05)/100)
print(sum(FPCA.out<0.05)/100)
print(sum(FdCor.out<0.05)/100)