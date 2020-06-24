# Cor   <- corr.fn(fields::rdist(S), theta = theta.C)
# eig.Sigma <- eigen(Cor, symmetric=TRUE)  # time consuming
# V <- eig.Sigma$vectors
# d <- eig.Sigma$values

# Dependencies: fields, truncnorm, mvtnorm


################################################################################
################################################################################
## Update Yi using full conditionals
## Y|β, range ~ N(Xβ,Σ), where Σ is the covariance matrix defined by 'range'.
## Zi|Y = Zi|Yi = I(Yi>0).
##
## When we update Yi, we immediately use its new value for other variables Yj
##
## i ................................. the location on which we are updating
## cen ................................. a (n.s x 1) vector of data for ONE day
## X ................................. a (n.s x 1) vector of current latent process
## X.s ................................. a (n.s x 1) vector of current latent process
## tau_sqd .............................. X*Beta, the regression parameters
## thresh.X

X.i.update<-function(i, cen, X, X.s, tau_sqd, thresh.X){
 
  ## Draw from full conditional using truncnorm
  if(!cen[i]) {a=thresh.X; b=Inf} else {a=-Inf; b=thresh.X}
  return(truncnorm::rtruncnorm(1, a=a, b=b, mean=X.s[i], sd=sqrt(tau_sqd)))
}
##
################################################################################
################################################################################


################################################################################
################################################################################
## The log likelihood of X.s, when it is conditioned on scaling factor R and Σ(λ,γ)
##
## X.s ............................... X, but without the measurement error
## R 
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##

X.s.likelihood.conditional<-function(X.s, R, V, d){
  if(any(X.s<R)) return(-Inf) else{
    X.s.to.Z <- qnorm(1-R/X.s)
    loglik <- -0.5*as.vector(eig2inv.quadform.vector(V, 1/d, X.s.to.Z))+0.5*sum(X.s.to.Z^2)-2*sum(log(X.s))-0.5*sum(log(d))+length(X.s)*log(R)
    return(loglik)
  }
}
##
################################################################################
################################################################################

################################################################################
##  Updates the X.s, for the generic Metropolis sampler
## Samples from the scaled Gaussian process (update the smooth process).
## The mixing distribution comes from from the Huser-wadsworth scale mixing distribution.
## The PROPOSAL will be the conjugate update from the model that treats the
## X process (i.e. X.s, but with measurement error) as the response, ignoring
## the marginal transformation.  Then the Metropolis part either rejects or
## accepts the entire draw.
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters (theta in dhuser.thibaud):
##                                     theta[1] = gamma (from H-O-T (2017))
##                                     theta[2] = beta (from H-O-T(2017))
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## theta.mix
## tau_sqd
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##

X.s.update.mixture.me.update.par.once.without.X.par <- function(R, X, X.s, cen,
                                                                prob.below, tau_sqd, V, d,
                                                                thresh.X=NULL){
  library(doParallel)
  library(foreach)
  n.s <- nrow(cen)
  n.t <- ncol(cen)
  registerDoParallel(cores=n.t)
  
  accepted <- matrix(0,nrow=n.s, ncol=n.t) 
  
  if(is.null(thresh.X)) thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
  
  Res<-foreach(t = 1:n.t, .combine = "cbind")%dopar%{
    res<- update_X_s_onetime(X = X[,t], X_s = X.s[,t], cen = cen[,t], prob_below = prob.below, 
                             delta = delta, tau_sqd = tau_sqd, thresh_X = thresh.X, R = R[t], V = V, d = d)
    c(res$X.s, res$accept)
  }
  
  X.s <- Res[1:n.s, ]
  accepted <-Res[(n.s+1):(2*n.s),]
  # O<-order(Res[2*n.s+1,])
  # X.s <- X.s[,O]
  # accepted <- accepted[,O]
  
  return(list(X.s=X.s, accepted=accepted))
}


#                                                                              #
################################################################################





################################################################################
## For the generic Metropolis sampler
## Samples from the scaling factors, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. R[t]
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error (vector)
## V, d
## delta
##
Rt.update.mixture.me.likelihood <- function(data, params, X.s, delta, 
                                            V=NULL, d=NULL) {
  R <- data
  Rt <- params
  if(Rt < 1)  return(-Inf)
  
  ll <- X.s.likelihood.conditional(X.s, Rt, V, d)
  return(as.vector(ll))
}

# Rt.update.mixture.me.likelihood(R, R[1], X.s[,1], delta, V, d)

#                                                                              #
################################################################################


################################################################################
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. delta (from Huser and Wadsworth 2017)
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## tau_sqd
##

delta.update.mixture.me.likelihood <- function(data, params, X, cen, 
                                               prob.below, tau_sqd) {
  R <- data
  delta <- params
  if(delta < 0 || delta > 1) return(-Inf)
  thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
  if(min(X[!cen])<thresh.X || max(X[cen])>thresh.X) return(-Inf)
  
  ll <- dhuser.wadsworth(R, delta, log=TRUE)
  
  return(ll)
}
# delta.update.mixture.me.likelihood(R, delta, Y, X.s, cen, prob.below, theta.gpd, tau)

#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the measurement error variance (on the X scale), for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
## Just a wrapper for marg.transform.data.mixture.me.likelihood
##
##   *********** If we do end up updating the prob.below parameter, this
##   *********** is a good place to do it.
##
## data............................... a n.t vector of scaling factors
## params............................. tau_sqd
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## delta
##
tau.update.mixture.me.likelihood <- function(data, params, X, X.s, cen, 
                                             prob.below, delta) {
  
  R <- data
  tau_sqd <- params
  
  thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
  if(min(X[!cen])<thresh.X || max(X[cen])>thresh.X) return(-Inf)
  
  ll <- sum(dnorm(X,mean=X.s,sd=sqrt(tau_sqd),log = 1))
  return(ll)
}
#tau.update.mixture.me.likelihood(R, tau, Y, X.s, cen, prob.below, delta, theta.gpd)

#                                                                              #
################################################################################




################################################################################
## Update covariance parameters. For the generic Metropolis sampler
## Samples from the parameters of the underlying Gaussian process, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters:
##                                     lambda = params[1]
##                                     gamma  = params[2]
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## R
## S
## V, d
##
theta.c.update.mixture.me.likelihood <- function(data, params, X.s, S, 
                                                 V=NULL, d=NULL) {
  
  library(doParallel)
  library(foreach)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  n.t <- ncol(X.s)
  registerDoParallel(cores=n.t)
  
  R <- data
  range <- params[1]
  nu <- params[2]
  # if(lambda<0 || gamma<0 || gamma>2)  return(-Inf)
  
  if(is.null(V)){
    Cor   <- corr.fn(rdist(S), theta=c(range,nu))
    eig.Sigma <- eigen(Cor, symmetric=TRUE)
    V <- eig.Sigma$vectors
    d <- eig.Sigma$values
  }
  
  ll<-foreach(i = 1:n.t, .combine = "c") %dopar% {
    X.s.likelihood.conditional(X.s[,i], R[i], V, d)
    # dmvn.eig(qnorm(1-R[i]/X.s[, i]), V = V, d.inv = 1/d)
  }
  
  return(sum(ll))
}

range.update.mixture.me.likelihood <- function(data, params, X.s, S, 
                                               nu, V=NULL, d=NULL) {
  
  library(doParallel)
  library(foreach)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  n.t <- ncol(X.s)
  registerDoParallel(cores=n.t)
  
  R <- data
  range <- params[1]
  # nu <- params[2]
  # if(lambda<0 || gamma<0 || gamma>2)  return(-Inf)
  
  if(is.null(V)){
    Cor   <- corr.fn(rdist(S), theta=c(range,nu))
    eig.Sigma <- eigen(Cor, symmetric=TRUE)
    V <- eig.Sigma$vectors
    d <- eig.Sigma$values
  }
  
  ll<-foreach(i = 1:n.t, .combine = "c") %dopar% {
    X.s.likelihood.conditional(X.s[,i], R[i], V, d)
    # dmvn.eig(qnorm(1-R[i]/X.s[, i]), V = V, d.inv = 1/d)
  }
  
  return(sum(ll))
}


#                                                                              #
################################################################################






################################################################################
## For the generic Metropolis sampler
## The likelihood for the range papameter (the smoothness nu is fixed)
##
## data............................... a (n.s x 1) vector of current latent process
## params............................. the range parameter
## X ................................. the design matrix,
##                                     sometimes equal to location matrix
## S ................................. location matrix
## nu ................................ the smoothness/roughness parameter
## Beta .............................. X*Beta, the regression parameters
##
range.update.likelihood <- function(data, params, X, S, nu, Beta, V=NULL, d=NULL) {
  
  Y <- data
  range <- params

  if(is.null(V)){
    Cor   <- corr.fn(fields::rdist(S), theta = c(range, nu))
    eig.Sigma <- eigen(Cor, symmetric=TRUE)
    V <- eig.Sigma$vectors
    d <- eig.Sigma$values
  }
  
  Mean <- X%*%Beta
  ll <- dmvn.eig(R=Y, mean=Mean, V=V, d.inv=1/d)
  
  return(sum(ll))
}

##
################################################################################

