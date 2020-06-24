source("~/Desktop/Research/scalemixture_clipped/scalemix_utils.R")
source("~/Desktop/Research/scalemixture_clipped/scalemix_likelihoods.R")
source("~/Desktop/Research/scalemixture_clipped/scalemix_priors.R")
source("~/Desktop/Research/scalemixture_clipped/generic_samplers.R")
source("~/Desktop/Research/scalemixture_clipped/scalemix_sampler_02.R")

library(fields)   # For rdist

# ------------ 1. Simulation settings -------------
n.s <- 100        # Number of sites
n.t <- 40         # Number of time points
tau <- 4          # Nugget SD
delta <- 0.53      # For R
range <- 1        # Matern range
nu <-  3/2        # Matern smoothness



# -------------- 2. Generate fake data -----------------
S     <- cbind(runif(n.s, 0, 5), runif(n.s, 0, 5))
# Cor   <- corr.fn(rdist(S), lambda = lambda, gamma = gamma)
Cor   <- corr.fn(rdist(S), theta=c(range,nu))
eig.Sigma <- eigen(Cor, symmetric=TRUE)
V <- eig.Sigma$vectors
d <- eig.Sigma$values

# set.seed(3333)
u<-rep(NA,n.t)
R<-rep(NA,n.t)
for(t in 1:n.t){
  u[t] <-runif(1,0,1)
  R[t]<-pow(1/(1-u[t]),delta/(1-delta))
}


X <- matrix(NA, n.s, n.t)
X.s <- matrix(NA, n.s, n.t)
for(t in 1:n.t) {
  Z.t<-eig2inv.times.vector(V, sqrt(d), rnorm(n.s))
  Z.to.W.s<-1/(1-pnorm(Z.t))
  X.s[ ,t] <- R[t]*Z.to.W.s
  X[ ,t] <- X.s[ ,t] + sqrt(tau)*rnorm(n.s)
}


prob.below <- 0.95

thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau, delta = delta) 
sum(X < thresh.X) / length(X)
cen <- X < thresh.X  



## --------------- 3. Saving initial values -------------------
initial.values <- list(delta = delta-0.02, tau=tau+1, theta.c=c(range,nu), prob.below=prob.below, X.s=X.s, R=R, S=S)
n.updates <- 50000
thin <- 10
echo.interval <- 50
true.params <- list(delta = delta, tau=tau, theta.c=c(range,nu), prob.below=prob.below, X.s=X.s, R=R, S=S)
save(X, thresh.X, cen, true.params, file="Initial.RData")




## --------------- 4. Running Metropolis -------------------
Res <- scalemix.sampler.02(S=S, cen=cen,
                    initial.values=initial.values,
                    n.updates=n.updates, thin=thin,
                    experiment.name="Huser-wadsworth-clipped",
                    echo.interval=echo.interval,
                    sigma.m=NULL, prop.Sigma=NULL,
                    true.params=true.params, sd.ratio=NULL, lower.prob.lim=0.5)

