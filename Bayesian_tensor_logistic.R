PGtensor <- function(Z,X,Ylabel, n=rep(1, length(Ylabel)),rank = 3,nsweep = 1e3, nskip = 3, a.lam, b.lam, phi.alpha, scale = TRUE) {
  library(base)
  library(statmod)
  library(MASS)
  library(matlib)
  library(stats)
  library(plyr)
  library(GIGrvg)
  library(gtools)
  library(coda)
  library(fields)
  library(dplyr)
  library(BayesLogit)
  
  ## data input
  # N : total number of subjects
  # P: Size of tensor images
  # D: Dimension of tensor images (2D/3D)
  # Z : N * pgamma matrix of demographic covariates
  # X  : N *  P[1] * P[2] tensor images 
  # Ylabel : N binary response
  # n_i = 1
  
  Yi = Ylabel
  N = length(Yi)
  P = dim(X)[-1]
  D = length(dim(X)) - 1
  pgamma <- ncol(Z) 
  
  # #### standarize ####
  if (is.null(Z)) {
    mz = NA; sz <- NA} else {
      mz = colMeans(Z); sz <- rep(1,pgamma)
    }
  if (!is.null(Z)) {
    if (!is.matrix(Z)) {Z <- as.matrix(Z)}
  }
  Zt <- Z
  obs <- as.numeric(Yi); sy=1; my=0

  if(scale){     ##centering & scale to have unity range
    if (!is.null(Z)) {
      mz <- colMeans(Z)
      sz <- apply(Z, 2, function(z) diff(range(z))); sz[sz==0] = 1
      Zt <- Z
      for (jj in 1:pgamma) {
        Zt[,jj] <- (Z[,jj] - mz[jj]) / sz[jj]
      }
    }
    
    Xt <- 0 * X
    mx <- apply(X, c(2:(D+1)), function(z) mean(z,na.rm=T))
    sx <- apply(X, c(2:(D+1)), function(z) diff(range(z,na.rm=T)))
    sx[sx == 0] <- 1
    if (D==2) {
      for(jj in 1:nrow(X)) Xt[jj,,] <- (X[jj,,] - mx) / sx
    } else if (D==3) {
      for(jj in 1:nrow(X)) Xt[jj,,,] <- (X[jj,,,] - mx) / sx
    }
  } else{
    ## do nothing;
    if (!is.null(Z)) {
      mz <- rep(0, pgamma)
      sz <- rep(1, pgamma)
      Zt <- Z
    }
    mx <- array(0, dim = dim(X)[-1])
    sx <- array(1, dim = dim(X)[-1])
    Xt <- X
  }
  
  #### MCMC setup ####
  
  #require(glmnet)
  x.train.nona = Xt
  x.train.nona[is.na(Xt)] = 0
  
  ## hyper-par initialize 
  if(missing(a.lam)) a.lam <- rep(3, rank)
  if(missing(b.lam)) b.lam <- (a.lam)**(1/(2 * D))
  if(missing(phi.alpha)) phi.alpha <- rep(1 / rank, rank)
  
  phi.a0 <- sum(phi.alpha)
  a.vphi <- phi.a0
  b.vphi <- phi.alpha[1] * rank^(1/D)
  
  s0 = 1; a.t = 0.1; b.t = 1 #b.t = 2.5/2 * s0^2
  #s0 = 1; a.t = 2.5/2; b.t = 2.5/2 * s0^2
  
  phi   <- rdirichlet(1, phi.alpha)
  varphi <- rgamma(1, a.vphi, b.vphi)
  tau.r <- phi * varphi
  
  ## Storage/Posterior Quantities
  if (is.null(Z)) {
    gam <- 0
  } else {
    gam = rep(0, pgamma) 
  }
  
  lambda <- matrix(rgamma(rank*D, a.lam[1], b.lam[1]), rank, D)
  omega <- lapply(1:D, function(x) array(rexp(rank*P[x],.5*(a.lam[1]/b.lam[1])),dim=c(rank,P[x])))
  beta <- lapply(1:D, function(x) array(rnorm(rank*P[x]),dim=c(rank,P[x]))) #D*rank*P_d
  
  ### Polya-gamma Parameters
  # Yi is 1/0 response, n is param from binom dist
  kappa = (Yi-1/2)*n
  w = rep(0,N)
  
  tens.mean = getmean(x.train.nona, beta, rank) # <X,B>
  if (!is.null(Z)) {
    pred.mean <- Zt %*% gam
  } else {
    pred.mean <- as.matrix(rep(0,length(Yi)))
  }
  yest = pred.mean+tens.mean 
  
  # Storage/Posterior Quantities
  w.store = matrix(0,N,nsweep)
  mu_store = matrix(0,N,nsweep)
  alpha.store <- rep(NA,nsweep)
  gam.store   <- array(data=NA,dim=c(nsweep,pgamma))
  phi.store   <- array(data=NA,dim=c(nsweep,rank))
  varphi.store <- array(data=NA,dim=c(nsweep,1))
  beta.store  <- lapply(1:nsweep, function(x) lapply(1:D, function(y) array(dim=c(rank,P[y]))))
  omega.store <- lapply(1:nsweep, function(x) lapply(1:D, function(y) array(dim=c(rank,P[y]))))
  lambda.store   <- array(data=NA,dim=c(nsweep,rank,D))
  hyppar.store <-array(data=NA,dim=c(nsweep,rank,2))
  
  # create a data frame
  par.grid <- expand.grid(alam=seq(2.1,D+1,length.out=5), 
                          zeta=seq(0.5,ceiling(10 * rank**(1/(2*D))/2)/10,length.out=5))
  #alpha.grid <- exp(log(rank) * seq(-d, -1/(2*d), length.out = 10)); M <- 10
  #alpha.grid <- exp(log(rank) * seq(-d, -0.1, length.out = 10)); M <- 20
  alpha.grid <- seq(rank**(-D), rank**(-0.1), length.out = 10); M <- 20
  score.store <- array(data=NA, dim=c(nsweep,length(alpha.grid)))
  
  ###################################### MCMC #########################################
  tt <- Sys.time()
  for (sweep in 1:nsweep) {
    
    ##################################################################
    
    # Sample w, the Polya-gamma parameter
    psi = drop(yest)
    for (sub in 1:N) {
      w = rpg.devroye(N, n, psi)
    }
    w.store[,sweep] = w
    
    ## Sample Gamma
    if (!is.null(Z)) {
      Z_w <- matrix(NA, N, pgamma) 
      for (i in 1:N) {
        Z_w[i, ] = Zt[i,]*sqrt(w[i]) #G
      }
      ZZ_w = crossprod(Z_w, Z_w) # G^T*G
      
      Sig.g <- chol2inv(chol(diag(pgamma) + ZZ_w))
      mu.g <- Sig.g %*% (crossprod(Zt*w, (kappa/w-tens.mean))) # incorporate rho in mean
      gam <- mu.g + chol(Sig.g) %*% rnorm(pgamma)
    } else {
      gam <- 0
    }
    ## Update pred.mean
    if (!is.null(Z)) {
      pred.mean <- Zt %*% gam
    } else {
      pred.mean <- as.matrix(rep(0,length(Yi)))
    }
    
    ## update (a.lam, b.lam)
    Cjr <- sapply(1:rank, function(rr){bb <- sapply(1:D, function(jj) sum(abs(beta[[jj]][rr,]))); bb <- bb / sqrt(tau.r[rr]); return(bb)})
    mfun <- function(z,rank){o <- sapply(1:D, function(x) return(lgamma(z[1]+P[x]) - lgamma(z[1]) + z[1]*log(z[2]*z[1]) - (z[1]+P[x])*log(z[2]*z[1]+Cjr[x,rank]))); return(sum(o))}
    ll <- sapply(1:rank, function(rr) apply(par.grid, 1, mfun, rank = rr))
    par.wt <- apply(ll, 2, function(z) return(exp(z - logsum(z))))
    ixx <- apply(par.wt, 2, sample, x = c(1:nrow(par.grid)), size = 1, replace = F)
    for(rr in 1:rank){
      a.lam[rr] <- par.grid[ixx[rr],1]
      b.lam[rr] <- par.grid[ixx[rr],2] * a.lam[rr]
    }
    
    ## update (alpha, phi, varphi)
    draw.phi_tau <- function(alpha){
      len <- length(alpha)
      
      m.phialpha <- rep(alpha[1], rank)
      m.phia0 <- sum(m.phialpha)
      m.avphi <- m.phia0
      ## assumes b.vphi const (use: alpha 1 / R)
      
      Cr <- sapply(1:rank, function(rr) {bb <- sapply(1:D, function(jj) crossprod(beta[[jj]][rr,], 
                                                                                  diag(1/omega[[jj]][rr,]) %*% beta[[jj]][rr,])); return(bb)})
      score.fn <- function(phi.alpha, phi.s, varphi.s, Cstat){
        ldirdens <- function(v, a){
          c1 <- lgamma(sum(a))
          c2 <- sum(lgamma(a))
          return((c1-c2) + sum((a-1) * log(v)))
        }
        ldir <- apply(phi.s, 1, ldirdens, a = phi.alpha)
        
        lvarphi <- dgamma(varphi.s, sum(phi.alpha), b.vphi, log = T)
        
        dnorm.log <- -rowSums(Cstat) / (2 * varphi.s) -(sum(P)/2) * sapply(1:length(varphi.s), function(ii)return(sum(log(varphi.s[ii] * phi.s[ii,]))))
        return(dnorm.log + ldir + lvarphi)
      }
      
      phi <- NULL; varphi <- NULL; scores <- NULL
      if(len > 1){
        phi <- matrix(0, M*length(alpha.grid), rank)
        varphi <- matrix(0, M*length(alpha.grid), 1)
        Cstat <- matrix(0, M*length(alpha.grid), rank)
        scores <- list()
        
        ## get reference set
        for(jj in 1:len){
          m.phialpha <- rep(alpha[jj], rank)
          m.phia0 <- sum(m.phialpha)
          m.avphi <- m.phia0
          
          ## draw phi
          Cr1 <- colSums(Cr)
          phi.a <- sapply(1:rank, function(rr){rgig(M,m.phialpha[rr]-sum(P)/2,Cr1[rr],2*b.vphi)})
          phi.a <- t(apply(phi.a, 1, function(z)return(z / sum(z)))) ## [M x rank]
          
          ## draw varphi ##colSums(Cr / t(replicate(d, z)))
          Cr2 <- t(apply(phi.a, 1, function(z)return(Cr1 / z)))
          varphi.a <- apply(Cr2, 1, function(z)return(rgig(1, m.avphi-rank*sum(P)/2, sum(z), 2*b.vphi)))
          phi[seq((jj-1)*M+1, jj*M), ] <- phi.a
          varphi[seq((jj-1)*M+1, jj*M)] <- varphi.a
          Cstat[seq((jj-1)*M+1, jj*M), ] <- Cr2                
        }
        scores <- lapply(alpha.grid, function(z)return(score.fn(rep(z,rank), phi, varphi, Cstat)))
        lmax <- max(unlist(scores))
        scores <- sapply(scores, function(z) return(mean(exp(z - lmax))))
      } else{
        ## draw phi
        Cr1 <- colSums(Cr)
        phi <- sapply(1:rank, function(rr){rgig(1,m.phialpha[rr]-sum(P)/2,Cr1[rr],2*b.vphi)})
        phi <- phi / sum(phi)
        
        ## draw varphi
        Cr2 <- Cr1 / phi
        varphi <- rgig(1, m.avphi-rank*sum(P)/2, sum(Cr2), 2*b.vphi)
      }
      return(list(phi=phi, varphi=varphi, scores=scores))
    }
    
    ## sample astar
    o <- draw.phi_tau(alpha.grid)
    astar <- sample(alpha.grid, size = 1, prob = o$scores)
    score <- o$scores/sum(o$scores)
    #cat(sprintf('scores: %s\n', paste(round(score <- o$scores/sum(o$scores),2),collapse = ', ')))
    score.store[sweep,] <- score
    
    ## sample (phi, varphi)
    o <- draw.phi_tau(astar)
    phi <- o$phi; varphi <- o$varphi
    tau.r <- varphi * phi
    phi.alpha <- rep(astar, rank); phi.a0 <- sum(phi.alpha); a.vphi <- phi.a0
    
    ###################### update rank specific params ########################
    for(r in 1:rank){
      for(j in 1:D){
        tens.mu.r <- getmean(x.train.nona, beta, rank, r)
        betj <- getouter_list(lapply(beta[-j],function(x) x[r,]))
        H <- matrix(NA, N, P[j]) # H
        H_w <- matrix(NA, N, P[j])
        
        for(i in 1:N) {
          if (D==2) {
            H[i, ] <- apply(x.train.nona[i,,], j, function(x) return(sum(x * betj)))
            H_w[i,] <- H[i, ]*sqrt(w[i])
          } else if (D==3) {
            H[i, ] <- apply(x.train.nona[i,,,], j, function(x) return(sum(x*betj)))
            H_w[i,] <- H[i, ]*sqrt(w[i])
          }
        }
        HH_w <- crossprod(H_w,H_w) #HH_w: H*Omega*H
        
        # posterior covariance
        K <- chol2inv(chol(HH_w + diag(1/omega[[j]][r,]) / tau.r[r]))
        
        # mm: y tilda
        mm = kappa/w - pred.mean - tens.mu.r 
        
        # posterior mean
        bet.mu.jr <- K %*% crossprod(H, mm*w)
        
        # update beta_jr
        beta[[j]][r,] <- bet.mu.jr + chol(K) %*% rnorm(P[j])
        
        ## update lambda.jr
        lambda[r,j] <- rgamma(1, a.lam[r] + P[j], b.lam[r] + sum(abs(beta[[j]][r,])) / sqrt(tau.r[r]))
        
        ## update omega.jr
        omega[[j]][r,] <- sapply(1:P[j], function(kk) rgig(1, 1/2,beta[[j]][r,kk]^2 / tau.r[r], lambda[r,j]^2))
        #omega[r,j,] <- sapply(1:p, function(kk){a <- lambda[r,j]^2; b <- beta[[r]][kk,j]^2 / tau.r[r]; map <- besselK(sqrt(a*b),0.5 + 1) / besselK(sqrt(a*b), 0.5) * sqrt(b / a); return(map)})
      }
    }
    
    tens.mean = getmean(x.train.nona, beta, rank)
    if (!is.null(Z)) {
      pred.mean <- Zt %*% gam
    }
    
    cat(sprintf('sweep: %s\n', sweep))
    
    ## store params
    beta.store[[sweep]] <- beta
    if (!is.null(Z)) {
      gam.store[sweep,] <- gam
    } else {
      gam.store[sweep] <- gam
    }
    alpha.store[sweep] <- astar ## not intercept
    phi.store[sweep,] <- phi
    varphi.store[sweep,] <- varphi
    omega.store[[sweep]] <- omega
    lambda.store[sweep,,] <- lambda
    sapply(1:rank, function(rr) hyppar.store[sweep,rr,] <<- c(a.lam[rr], b.lam[rr]))
    }
  
  
  tt <- abs(tt - Sys.time())
  cat('Time out: ', tt, '\n')
  
  #### plotting omited####
  
  #### finalize ####    
  out <- list(nsweep = nsweep, 
              rank = rank, 
              mu_store = mu_store,
              P = P, 
              D = D, 
              w.store = w.store,
              my = my, sy = sy, mz = mz, sz = sz, mx = mx, sx = sx, gam.store = gam.store, beta.store = beta.store, time = tt)
  
  class(out) <- "tensor.reg"
  return(out)
  
}


#### aux functions ####

getouter_list <- function(bet) {
  D <- length(bet)
  if (D==1) {
    return(bet[[1]])
  }
  if (D==2) {
    return(outer(bet[[1]], bet[[2]]))
  }
  else {
    return(outer(getouter_list(bet[1:(D-1)]), bet[[D]]))
  }
}

TP.rankR <- function(X.allr) { 
  R <- ncol(X.allr[[1]])
  if (is.null(R)) {
    return(getouter_list(X.allr))
  } else {
    Y <- array(0, dim=c(as.numeric(lapply(X.allr, function(x) length(x[,1])))))
    for (r in c(1:R)) {
      Y <- Y + getouter_list(lapply(X.allr, function(x) x[,r]))
    }
    return(Y)
  }
}

getmean <- function(X,beta,rank,rank.exclude=NULL){    
  idx <- setdiff(1:rank,rank.exclude)
  B <- Reduce('+', lapply(idx,function(r) getouter_list(lapply(beta,function(x) x[r,]))))
  mu.B <- apply(X, 1, function(xx, bb) sum(xx * bb), bb = B)
  return(mu.B)
}

tensor.mean = function(x,n) {
  Reduce("+", x)/n
}

logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))

getBeta_mcmc <- function(beta.store) {
  nsweep <- length(beta.store)
  D <- length(beta.store[[1]])
  rank <- nrow(beta.store[[1]][[1]])
  P <- sapply(1:D,function(x) ncol(beta.store[[1]][[x]]))
  Beta_mcmc <- array(dim=c(nsweep,prod(P)))
  for (i in 1:nsweep) {
    coef <- rep(0,prod(P))
    for (r in 1:rank) {
      coef <- coef + c(getouter_list(lapply(beta.store[[i]],function(x) x[r,])))
    }
    Beta_mcmc[i,] <- coef
  }
  quantile(Beta_mcmc)
  return(Beta_mcmc)
}

add <- function(x) Reduce("+", x)

uncollapse <- function(str, collapse = "", mode = "character"){
  a <- unlist(strsplit(str, collapse))
  mode(a) <- mode
  return(a)
}

rmse <- function(y, yhat) {
  return(sqrt(mean((y-yhat)^2,na.rm=T)))
}

logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))


# #### sample on simulated data ####
# library(fields)
# set.seed(123)
# # Generate tensor signal
# Beta_tens = matrix(0, 48, 48)
# for (i in 15:40) {
#   for (j in 10:35) {
#     Beta_tens[i,j] = 1
#   }
# }
# # Generate 2D images and binary response
# N = 500; p = c(48,48); rank = 3
# X = array(rnorm(N*prod(p)),dim=c(N,p)) # simulated image
# Y = sapply(1:N,function(x) sum(X[x,,]*Beta_tens,na.rm=T))
# hist(Y)
# p = 1/(1+exp(-Y))
# Ylabel <- rbinom(n = N, size = 1, prob = p)
# 
# train_index = sort(sample(1:N, 0.7*N, replace = FALSE))
# x.train = X[train_index,,]
# y.train = Ylabel[train_index]
# 
# burnin = 1000; nsweep = 3000
# sim = PGtensor(Z = NULL, x.train, y.train,n=rep(1, length(y.train)), nsweep=nsweep, rank=3, scale=T)
# 
# # tensor_est: estimated tensor coefficient
# tensor = getBeta_mcmc(sim$beta.store);tensor_est = apply(tensor[burnin:nsweep,],2,mean)*sim$sy/sim$sx
# rmse(c(Beta_tens),c(tensor_est))
# cor(c(Beta_tens),c(tensor_est))
# # plot estimated tensor coefficient
# image.plot(tensor_est, col = gray.colors(25, start =1, end = 0), axes = F)
# mtext(text=seq(10,50,10), side=2, line=0.3, at=seq(10,50,10)/48, las=1, cex=0.8)
# mtext(text=seq(10,50,10), side=1, line=0.3, at=seq(10,50,10)/48, las=2, cex=0.8)
# # get misclassification rate
# x.test = X[-train_index,,]
# y.test = Ylabel[-train_index]
# post.tens.mean.test = array(0, N-length(train_index))
# for (j in 1:(N-length(train_index))) {post.tens.mean.test[j] = c(tensor_est)%*%c(x.test[j,,])}
# post.mui.test = matrix(post.tens.mean.test, ncol = 1)
# clust.test = rep(0, N-length(train_index))
# clust.test[post.mui.test>0] = 1
# clust.test[post.mui.test<=0] = 0
# missclassrate=1-sum(clust.test == y.test)/length(clust.test)
# # get f1 score
# TP = 0; FP = 0; FN = 0
# for (j in 1:length(clust.test)) {
#   if (clust.test[j] == 1 && y.test[j] == 1) {TP = TP+1}
#   else if (clust.test[j] == 1 && y.test[j] == 0) {FP = FP+1} 
#   else if (clust.test[j] == 0 && y.test[j] == 1) {FN = FN+1}
# }
# f1score = TP/(TP+(FP+FN)/2)


