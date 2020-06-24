library(fields)
library(scales)

###################################################################################
## Main sampler

## Y ........................................... censored observations on GPD scale
## S ............................................................. location indices
## cen ........................................................... indicator matrix
## thresh .................................................... threshold on Y scale
## initial.values .................. a list: delta, rho, tau, theta.gpd, prob.below
##                                             X.s, R
## n.updates .................................................... number of updates
## thin ............................................. number of runs in each update
## experiment.name
## echo.interval ......................... echo process every echo.interval updates
## sigma.m
## prop.Sigma
## true.params ..................... a list: delta, rho, tau, theta.gpd, prob.below
##                                              X.s, R
## lower.prob.lim ......................................... for updating prob.below

scalemix.sampler.02 <- function(S, cen, initial.values,
                                n.updates, thin=10,
                                experiment.name="Huser-wadsworth",
                                echo.interval=50,
                                sigma.m=NULL, prop.Sigma=NULL, 
                                true.params=NULL, sd.ratio=NULL, lower.prob.lim=0.5) {
  
  #library(doParallel)
  #library(foreach)
  save.bit <- TRUE
  
  # Constants to control how many Metropolis steps get executed at each
  # iteration of the main loop
  # n.metr.updates.prob.below <- 4
  n.metr.updates.delta <- 4
  n.metr.updates.theta.c <- 4
  n.metr.updates.tau <- 4
  n.metr.updates.R <- 4
  
  # Constants to control adaptation of the Metropolis sampler
  c.0 <- 1
  c.1 <- 0.8
  k <- 3  # the iteration offset
  metr.opt.1d <- 0.41
  metr.opt.2d <- 0.35
  
  # Hyper parameters for the prior of the mixing distribution parameters and 
  # the correlation parameters
  # hyper.params.prob.below <- c(0, 10)
  hyper.params.delta <- c(0,1)
  hyper.params.tau <- 100 #c(0.1,0.1)
  hyper.params.theta.c <- c(1,1)
  hyper.params.range <- c(0.5,1.5) # in case where roughness is not updated
  
  # A small number
  eps <- 1e-06
  
  # Bookkeeping
  n.s <- nrow(cen)
  n.t <- ncol(cen)
  h <- rdist(S)
  diag(h)  <- 0
  #registerDoParallel(cores=n.t)
  
  # Load initial values
  delta <- initial.values$delta
  tau <- initial.values$tau
  theta.c <- initial.values$theta.c
  prob.below <- initial.values$prob.below
  X.s <- initial.values$X.s
  R <- initial.values$R
  
  X.s.accept <- rep(1, n.t)
  
  
  
  #**********   Eigendecomposition of the correlation matrix   ****************#
  Sigma <- corr.fn(h, theta=theta.c)
  eig.Sigma <- eigen(Sigma, symmetric=TRUE)
  V <- eig.Sigma$vectors
  d <- eig.Sigma$values
  thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, delta = delta)
  
  # Initialize trace objects
  accepted <- rep(0, n.updates)
  save.num <- n.t; subset.replicates <- 1:n.t
  if(n.t>60) {save.num <- 30;subset.replicates<-ceiling(seq(1,n.t,length.out = save.num))}
  # prob.below.trace <- rep(NA, n.updates)
  X.s.trace <- array(NA, dim=c(n.updates, n.s, save.num))
  X.trace <- array(NA, dim=c(n.updates, n.s, save.num))
  X.s.accept.trace <- matrix(NA, n.updates, n.t)
  tau.trace <- rep(NA, n.updates)
  delta.trace <- rep(NA, n.updates)
  theta.c.trace <- matrix(NA, n.updates, 2)
  R.trace <- matrix(NA, n.updates, n.t)
  
  
  # Column names for trace objects
  colnames(theta.c.trace)   <- c("range", "nu")
  
  # Fill in trace objects with initial values
  X.s.trace[1, , ] <- X.s[,subset.replicates]
  X.trace[1, , ] <- X[,subset.replicates]
  X.s.accept.trace[1, ] <- X.s.accept
  tau.trace[1] <- tau
  delta.trace[1] <- delta
  theta.c.trace[1, ] <- theta.c
  # prob.below.trace[1] <- prob.below
  R.trace[1, ] <- R
  
  
  # For tuning Metropolis updates of theta
  # if (is.null(sigma.m$prob.below)) sigma.m$prob.below <- 1
  # if (is.null(sigma.m$X.s)) sigma.m$X.s <- matrix(2.4^2,nrow=n.s,ncol=n.t)
  if (is.null(sigma.m$delta)) sigma.m$delta <- 2.4^2 #0.05^2
  if (is.null(sigma.m$theta.c)) sigma.m$theta.c<- (2.4/2)^2
  if (is.null(sigma.m$range)) sigma.m$range<- 2.4/2  # in case where roughness is not updated
  if (is.null(sigma.m$tau)) sigma.m$tau <- 2.4^2
  if (is.null(sigma.m$R)) sigma.m$R <- rep(2.4^2, n.t)
  
  if(is.null(prop.Sigma$theta.c))  prop.Sigma$theta.c<-diag(2)

  
  r.hat.delta <- NA
  r.hat.theta.c <- NA
  # r.hat.range <- NA  # in case where roughness is not updated
  r.hat.prob.below <- NA
  r.hat.tau <- NA
  # r.hat.R <- resp(NA, n.t)
  
  
  for (i in 2:n.updates) {
    
    ################################################################
    ## Update Metropolis adaptation parameters
    ################################################################
    gamma2 <- c.0 / (i + k)^(c.1)
    gamma1 <- 1 / (i + k)^(c.1)
    # Accepted<-matrix(0, nrow=n.s, ncol=n.t) #For X.s
    
    for (j in 1:thin) {
      
      ################################################################
      ## Update X -- branching: use X.s
      ################################################################
      library(doParallel)
      library(foreach)
      registerDoParallel(cores=n.t)
      X<-foreach(t = 1:n.t, .combine = "cbind") %dopar% {
        tmp <- X[,t]
        for(loca in 1:n.s){
          tmp[loca] <- X.i.update(loca, cen[,t], tmp, X.s[,t], tau, thresh.X)      
        }
        tmp      
      }
      
      
      ################################################################
      ## Update X.star  *Parallel
      ################################################################
      X.s.res<-X.s.update.mixture.me.update.par.once.without.X.par(R=R,X=X, X.s=X.s, cen=cen, prob.below=prob.below, 
                                                                   tau_sqd=tau, V=V, d=d, thresh.X=thresh.X)
      X.s <- X.s.res$X.s
      X.s.accept <- apply(X.s.res$accepted==1,2,any)
      # Accepted <- tryCatch(Accepted + X.s.res$accepted, error=function(e){cat("The dim of X.s.res$accepted is (",dim(X.s.res$accepted),")","\n")})
      
      
      
      ################################################################
      ## Update theta.c
      ################################################################
      # metr.out.theta.c <- static.metr(z = R, starting.theta = theta.c,
      #                   likelihood.fn = theta.c.update.mixture.me.likelihood, 
      #                   prior.fn = half.cauchy.half.cauchy, hyper.params = hyper.params.theta.c, 
      #                   n.updates = n.metr.updates.theta.c, prop.Sigma = prop.Sigma$theta.c, 
      #                   sigma.m=sigma.m$theta.c, verbose=FALSE,
      #                   X.s = X.s, S = S)
      # r.hat.theta.c <- metr.out.theta.c$acc.prob
      # theta.c <- metr.out.theta.c$trace[n.metr.updates.theta.c, ]
      # sigma.m$theta.c <- exp(log(sigma.m$theta.c) + gamma2*(r.hat.theta.c - metr.opt.2d))
      # prop.Sigma$theta.c <- prop.Sigma$theta.c + gamma1*(cov(metr.out.theta.c$trace)-prop.Sigma$theta.c)
      
      
      ## in case where roughness is not updated
      metr.out.theta.c <- static.metr(z = R, starting.theta = theta.c[1],
                                      likelihood.fn = range.update.mixture.me.likelihood,
                                      prior.fn = interval.unif, hyper.params = hyper.params.range,
                                      n.updates = n.metr.updates.theta.c, prop.Sigma = 1,
                                      sigma.m=sigma.m$range, verbose=FALSE,
                                      X.s = X.s, S = S, nu = theta.c[2])
      r.hat.theta.c <- metr.out.theta.c$acc.prob
      theta.c[1] <- metr.out.theta.c$trace[n.metr.updates.theta.c]
      sigma.m$range <- exp(log(sigma.m$range) + gamma2*(r.hat.theta.c - metr.opt.1d))
      
      
      ## Re-create covariance matrix and eigenvectors/eigenvalues
      if(r.hat.theta.c>0) {
        Sigma   <- corr.fn(h, theta=theta.c)
        eig.Sigma <- eigen(Sigma, symmetric=TRUE)
        V <- eig.Sigma$vectors
        d <- eig.Sigma$values
      }
      
      
      
      ################################################################
      ## Update R  *Parallel
      ################################################################
      Metr.R<-foreach(t = 1:n.t, .combine = "rbind") %dopar% {
        
        metr.out.R <- static.metr(z = R, starting.theta = R[t],
                                  likelihood.fn = Rt.update.mixture.me.likelihood, 
                                  prior.fn = huser.wadsworth.prior, hyper.params = delta,
                                  n.updates = n.metr.updates.R, prop.Sigma = 1, 
                                  sigma.m=sigma.m$R[t], verbose=FALSE,
                                  X.s = X.s[,t], delta = delta, V = V, d = d)
        c(metr.out.R$trace[n.metr.updates.R],
          exp(log(sigma.m$R[t]) + gamma2*(metr.out.R$acc.prob - metr.opt.1d)))
      }
      
      R <- Metr.R[, 1]
      sigma.m$R <- Metr.R[, 2]
      
      
      
      ################################################################
      ## Update tau
      ################################################################
      metr.out.tau <-  static.metr(z = R, starting.theta = tau, 
                                   likelihood.fn = tau.update.mixture.me.likelihood, 
                                   prior.fn = tau.sqd.prior, hyper.params = hyper.params.tau,
                                   n.updates = n.metr.updates.tau, prop.Sigma = 1, 
                                   sigma.m=sigma.m$tau, verbose=FALSE, 
                                   X = X, X.s = X.s, cen = cen, prob.below = prob.below, 
                                   delta = delta)
      tau <- metr.out.tau$trace[n.metr.updates.tau]
      r.hat.tau <- metr.out.tau$acc.prob
      sigma.m$tau <- exp(log(sigma.m$tau) + gamma2*(r.hat.tau - metr.opt.1d))
      
      
      
      ################################################################
      ## Update delta
      ################################################################
      metr.out.delta <- static.metr(z = R, starting.theta = delta, 
                                    likelihood.fn = delta.update.mixture.me.likelihood, 
                                    prior.fn = interval.unif, hyper.params = hyper.params.delta, 
                                    n.updates = n.metr.updates.delta, prop.Sigma = 1, 
                                    sigma.m=sigma.m$delta, verbose=FALSE, 
                                    X=X, cen = cen, prob.below = prob.below, tau_sqd = tau)
      delta <- metr.out.delta$trace[n.metr.updates.delta]
      r.hat.delta <- metr.out.delta$acc.prob
      if(sigma.m$delta>1e-4) sigma.m$delta <- exp(log(sigma.m$delta) + gamma2*(r.hat.delta - metr.opt.1d))
      
      
      
      ################################################################
      ## Update prob.below
      ################################################################
      #    metr.out.prob.below <- static.metr(R, logit.prob.below,
      #                                       logit.prob.below.update.mixture.me.likelihood,
      #                                       normal.scalar,
      #                                       hyper.params=hyper.params.prob.below,
      #                                       n.updates=n.metr.updates.prob.below,
      #                                       prop.Sigma=1,
      #                                       sigma.m=sigma.m$prob.below,
      #                                       Y=Y, X.s=X.s, cen=cen,
      #                                       theta.mix=theta.mix, theta.gaussian=theta.gaussian,
      #                                       theta.gpd=theta.gpd, lower.prob.lim=lower.prob.lim)
      #    logit.prob.below <- metr.out.prob.below$trace[n.metr.updates.prob.below]
      #    prob.below <- ilogit(logit.prob.below, c(lower.prob.lim, 1))
      #    r.hat.prob.below <- metr.out.prob.below$acc.prob
      #    sigma.m$prob.below <- exp(log(sigma.m$prob.below) +
      #                    gamma1*(r.hat.prob.below - metr.opt.1d))
      
      # Re-calculate the threshold on the X-scale    
      thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, delta = delta)
    }
    
    
    
    # # ---------------- Attempt to adapt proposal covariance ------------------
    # if (r.hat.theta.gpd > 0) {
    #   sd.ratio.hat <- sd(metr.out.theta.gpd$trace[ ,1]) / sd(metr.out.theta.gpd$trace[ ,2])
    # } else {
    #   sd.ratio.hat <- 1
    # }
    # sd.ratio <- exp(log(sd.ratio) + gamma1*(log(sd.ratio.hat) - log(sd.ratio)))
    # prop.Sigma$theta.gpd <-  matrix(c(1, prop.Sigma$gpd.corr/sd.ratio, prop.Sigma$gpd.corr/sd.ratio, 1/sd.ratio^2), 2, 2)
    # 
    # 
    # # ---------------- Adapt sigma.m$X.s ------------------
    # r.hat.X.s <- Accepted/thin
    # sigma.m$X.s <- exp(log(sigma.m$X.s) + gamma1*(r.hat.X.s - metr.opt.1d))
    
    
    # ------------------------- Fill in trace objects ---------------------------
    # prob.below.trace[i] <- prob.below
    X.s.trace[i, , ] <- X.s[,subset.replicates]
    X.trace[i, , ] <- X[,subset.replicates]
    X.s.accept.trace[i, ] <- X.s.accept
    theta.c.trace[i,]<-theta.c
    tau.trace[i] <- tau
    delta.trace[i] <- delta
    R.trace[i, ] <- R
    
    
    
    # ----------------------------- Echo progress --------------------------------
    cat("Done with", i, "updates,\n")
    if ((i %% echo.interval) == 0) {
      cat(" Acceptance rate for X^* is", mean(X.s.accept.trace[(i-echo.interval):i, ]), "\n")
      pdf(file=sprintf("%s_progress.pdf", experiment.name))
      par(mfrow=c(3,2))
      plot(theta.c.trace[,1], type="l", ylab=expression(range))
      if (!is.null(true.params)) abline(h=true.params$theta.c[1], lty=2, col=2, lwd=3)
      plot(theta.c.trace[,2], type="l", ylab=expression(nu))
      if (!is.null(true.params)) abline(h=true.params$theta.c[2], lty=2, col=2, lwd=3)
      plot(tau.trace, type="l", ylab=expression(tau^2))
      if (!is.null(true.params)) abline(h=true.params$tau, lty=2, col=2, lwd=3)
      plot(delta.trace, type="l", ylab=expression(delta))
      if (!is.null(true.params)) abline(h=true.params$delta, lty=2, col=2, lwd=3)
      # plot(prob.below.trace, type="l", ylab=expression(pi))
      # if (!is.null(true.params)) abline(h=true.params$prob.below, lty=2, col=2, lwd=3)
      par(mfrow=c(3, 3))
      for (j in 1:min(18, n.t)) {
        plot(R.trace[ , j], type="l", main=paste("Replication", j), 
             ylab=paste0("R[", j, "]"))
        if (!is.null(true.params)) abline(h=true.params$R[j], lty=2, lwd=3, col="gray80")
      }
      # For each year, plot the location with an exceedance that happens to be first in the matrix
      for (j in 1:min(18, save.num)) { 
        wh <- subset.replicates[j]
        plot.loc <- which(!cen[ ,wh])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylab=paste0("X^*[", plot.loc, ",", wh, "]"))
          if (!is.null(true.params)) abline(h=true.params$X.s[plot.loc, wh], lty=2, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      # For each year, plot the location with a non-exceedance that happens to be first in the matrix
      for (j in 1:min(18, save.num)) { 
        wh <- subset.replicates[j]
        plot.loc <- which(cen[ ,wh])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylim=c(min(X.s.trace[ , plot.loc, j], na.rm=TRUE), thresh.X),
               ylab=paste0("X^*[", plot.loc, ",", j, "]"))
          abline(h=thresh.X, lty=1, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      dev.off()
      
      # lines(curr.draw, col=alpha("gray80", 0.05))
      # # points(X.imp[ ,j], col=2, pch=20, cex=0.5)
      # Sys.sleep(0.1)
      state <- list(cen=cen, S=S, thresh=thresh,
                    experiment.name="Huser-wadsworth",
                    i=i, sigma.m=sigma.m, prop.Sigma=prop.Sigma, 
                    X=X, X.s=X.s, prob.below=prob.below,
                    delta=delta, R=R, tau=tau, theta.c=theta.c)
      out.obj <- list(cen=cen, S=S, thresh=thresh,
                      i=i,
                      X.s.trace=X.s.trace,
                      X.trace=X.trace,
                      X.s.accept.trace=X.s.accept.trace,
                      theta.c.trace=theta.c.trace,
                      tau.trace=tau.trace,
                      delta.trace=delta.trace,
                      # prob.below.trace=prob.below.trace,
                      R.trace=R.trace)
      save(state, out.obj, file=sprintf("%s_progress_%1.1s.RData", experiment.name, save.bit))
      save.bit <- !save.bit
    }
  }
  
  return(out.obj)
}
