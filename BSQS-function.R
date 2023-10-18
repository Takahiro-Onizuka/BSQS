library(MASS)
library(MCMCpack)
library(GIGrvg)
library(invgamma)
library(Bessel)
library(EigenR)

# -------------------------------------------------- #
#     Bayesian spatial quantile trend filtering      #
#               under horseshoe prior                #
# -------------------------------------------------- #

BSQTF.HS <- function(X, Y, D, pp=0.2, mc=1000, sparse=F, coordinatewise=F, n.print=F){
  
  if(n.print==F){TorF <- F}else{TorF <- T}
  
  # dataset
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- meany <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
    meany[i] <- mean(data$Y)
  }
  
  if(dim(D)[2] != length(location)){
    stop("The number of observation points and columns in D must be equal.")
  }
  
  
  ## settings
  m <- dim(D)[1]  # size of matrix
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  
  TD <- t(D)
  if(sparse){
    TD <- as(TD, "sparseMatrix")
    D <- as(D, "sparseMatrix")
  }
  
  
  ## MCMC box
  th.pos <- matrix(NA, mc, n)
  sig.pos <- c()
  rho.pos <- c()
  
  
  ## initial value
  th <- meany
  sig2 <- 1
  ww <- rep(1, m)   # local parameter in HS
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  tau2 <- 1   # global parameter in HS
  nu <- rep(1, m)   # latent parameter in HS
  xi <- 1   # latent parameter in HS
  vv <- rep(1, N)   # latent parameter in asymmetric Laplace,Exp(sig2)
  IGa <- 1   # paramater of sig2 prior
  IGb <- 1   # paramater of sig2 prior
  
  # acceptance rate
  ac.rho <- 0
  
  ## MCMC replications
  for(k in 1:mc){
    # theta
    if(coordinatewise==T){
      TDWD <- as.matrix(TD%*%diag(1/(ww*tau2))%*%D)
      A1 <- as.vector(sapply(split(vv^(-1), X), sum))/t2 + diag(TDWD)
      B1 <- as.vector(sapply(split((Y-psi*vv)/vv, X), sum))/t2
      for (i in 1:n) {
        ai <- A1[i]
        bi <- B1[i] - t(TDWD[i,-i])%*%th[-i]
        th[i] <- rnorm(1, bi/ai, sqrt(sig2/ai))
      }
      th.pos[k,] <- th
    }else{
      cc <- as.vector(sapply(split(vv^(-1), X), sum))/t2
      A <- solve(TD%*%diag(1/(ww*tau2))%*%D + diag(cc) )
      A <- (t(A) + A)/2
      B <- as.vector(sapply(split((Y-psi*vv)/vv, X), sum))/t2
      th <- mvrnorm(1, A%*%B, sig2*A)
      th.pos[k,] <- th
    }
    
    ## asymetric Laplace 
    # sig2
    eta <- as.vector(D%*%th)
    bb <- rep(0, N)
    for (i in 2:(n+1)) {
      bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-th[i-1]
    }
    cc <- sum(((bb-psi*vv)^2)/vv)/(2*t2) + sum(eta^2/ww)/(2*tau2) + sum(vv) + IGb
    sig2 <- rinvgamma(1, (n+3*N)/2+IGa, cc)   
    sig.pos[k] <- sig2
    
    # vv
    bb2 <- bb^2/(t2*sig2)
    cc <- 2/sig2 + psi^2/(t2*sig2)
    vv <- tapply(bb2, 1:N, rgig, n=1, lambda=1/2, psi=cc)
    vv[vv<lower] <- lower
    
    ## HS latent variables
    #ww
    bb <- eta^2/(2*tau2*sig2)+1/nu
    ww <- rinvgamma(m, 1, bb)
    ww[ww<lower] <- lower
    ww[ww>upper] <- upper
    
    #nu
    nu <- rinvgamma(m, 1/2, 1+1/ww)
    
    #tau2
    bbb <- sum(eta^2/ww)/(2*sig2) + 1/xi
    tau2 <- min(max(rinvgamma(1, (n+1)/2, bbb),10^(-5)), 10^5)   
    
    #xi
    xi <- rinvgamma(1, 1/2, 1+1/tau2)
    
    # print 
    if(TorF==T){ if(round(k/n.print)==(k/n.print)){ print(k) } }
  }
  
  return(th.pos)
}





# -------------------------------------------------- #
#     Bayesian spatial quantile trend filtering      #
#                under Laplace prior                 #
# -------------------------------------------------- #


BSQTF.Lap <- function(X, Y, D, pp=0.5, mc=1000, sparse=F, n.print=F){
  
  if(n.print==F){TorF <- F}else{TorF <- T}
  
  # dataset
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
  }
  
  if(dim(D)[2] != length(location)){
    stop("The number of observated locations and columns in D must be equal.")
  }
  
  ## settings
  m <- dim(D)[1]  # size of matrix
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  
  TD <- t(D)
  if(sparse){
    TD <- as(TD, "sparseMatrix")
    D <- as(D, "sparseMatrix")
  }
  
  
  ## MCMC box
  th.pos <- matrix(NA, mc, n)
  sig.pos <- c()
  
  
  ## initial value
  th <- y
  sig2 <- 1
  gam2 <- 1
  vv <- rep(1, N)   # latent parameter in asymmetric Laplace,Exp(sig2)
  ww <- rep(1, m)   #latent parameter in Laplace
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  d2 <- 1  # parameter of nu prior
  nu <- 1   # latent parameter in Laplace
  IGb <- 1   # paramater of sig2 prior
  IGa <- 1   # paramater of sig2 prior
  
  ## MCMC replications
  for(k in 1:mc){
    # theta
    cc <- as.vector(sapply(split(vv^(-1), X), sum))/t2
    A <- solve(TD%*%diag(1/ww)%*%D + diag(cc) )
    A <- (t(A) + A)/2
    B <- as.vector(sapply(split((Y-psi*vv)/vv, X), sum))/t2
    th <- mvrnorm(1, A%*%B, sig2*A)
    th.pos[k,] <- th
    
    ## asymetric Laplace 
    # sig2
    eta <- as.vector(D%*%th)
    bb <- rep(0, N)
    for (i in 2:(n+1)) {
      bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-th[i-1]
    }
    cc <- sum(((bb-psi*vv)^2)/vv)/(2*t2) + sum(eta^2/ww)/2 + sum(vv) + IGb
    sig2 <- rinvgamma(1, (n+3*N)/2+IGa, cc )
    sig.pos[k] <- sig2
    
    # vv
    bb2 <- bb^2/(t2*sig2)
    cc <- 2/sig2 + psi^2/(t2*sig2)
    vv <- tapply(bb2, 1:N, rgig, n=1, lambda=1/2, psi=cc)
    vv[vv<lower] <- lower
    
    ## Laplace latent variables
    #ww
    eta <- as.vector(D%*%th)
    bb3 <- eta^2/sig2
    ww <- tapply(bb3, 1:m, rgig, n=1, lambda=1/2, psi=gam2)
    ww[ww<lower] <- lower
    ww[ww>upper] <- upper
    
    #gam2
    gam2 <- rgig(1, lambda=m-1/2, chi=2/nu, psi=sum(ww))

    #nu
    nu <- rinvgamma(1, 1/2, 1/d2+1/gam2)
    
    # print 
    if(TorF==T){
      if(round(k/n.print)==(k/n.print)){ print(k) }
    }
  }

  return(th.pos)
}





# -------------------------------------------------- #
#         Bayesian spatial quantile smoothing        #
#                   under SAR prior                  #
# -------------------------------------------------- #

BSQS.SAR <- function(X, Y, W, pp=0.2, mc=1000, delta=0.2, sparse=F, coordinatewise=F, n.print=F){
  
  if(n.print==F){TorF <- F}else{TorF <- T}
  
  # dataset
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- meany <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
    meany[i] <- mean(data$Y)
  }
  
  if(dim(W)[2] != length(location)){
    stop("The number of observated lopcations and columns in W must be equal.")
  }
  
  
  ## settings
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  
  TW <- t(W)
  if(sparse){
    TW <- as(TW, "sparseMatrix")
    W <- as(W, "sparseMatrix")
  }
  
  # initial values
  rho <- 0.5
  IW <- as.matrix(diag(1,n)-rho*W)
  TIW <- as(t(IW),"sparseMatrix")
  IW <- as(IW,"sparseMatrix")
  iQ <- as.matrix(TIW%*%(IW))
  detiQ <- Eigen_det(asSparseMatrix(iQ))
  ac.rho <- 0
  
  ## MCMC box
  th.pos <- matrix(NA, mc, n)
  sig.pos <- c()
  rho.pos <- c()
  
  
  ## initial value
  th <- meany
  sig2 <- 1
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  tau2 <- 1   # 
  a_tau <- b_tau <- 0.1
  vv <- rep(1, N)   # latent parameter in asymmetric Laplace,Exp(sig2)
  IGa <- 1   # paramater of sig2 prior
  IGb <- 1   # paramater of sig2 prior
  
  
  ## MCMC replications
  for(k in 1:mc){
    # theta
    if(coordinatewise==T){
      TDWD <- iQ/tau2
      A1 <- as.vector(sapply(split(vv^(-1), X), sum))/t2 + diag(TDWD)
      B1 <- as.vector(sapply(split((Y-psi*vv)/vv, X), sum))/t2
      for (i in 1:n) {
        ai <- A1[i]
        bi <- B1[i] - t(TDWD[i,-i])%*%th[-i]
        th[i] <- rnorm(1, bi/ai, sqrt(sig2/ai))
      }
      th.pos[k,] <- th
    }else{
      cc <- as.vector(sapply(split(vv^(-1), X), sum))/t2
      A <- solve(iQ/tau2 + diag(cc) )
      A <- (t(A) + A)/2
      B <- as.vector(sapply(split((Y-psi*vv)/vv, X), sum))/t2
      th <- mvrnorm(1, A%*%B, sig2*A)
      th.pos[k,] <- th
    }
    
    ## asymetric Laplace 
      # sig2
      bb <- rep(0, N)
      for (i in 2:(n+1)) {
        bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-th[i-1]
      }
      cc <- sum(((bb-psi*vv)^2)/vv)/(2*t2) + sum(vv) + t(th)%*%iQ%*%th/tau2 + IGb
      sig2 <- rinvgamma(1, (n+3*N)/2+IGa, cc)   
      sig.pos[k] <- sig2
      
      # vv
      bb2 <- bb^2/(t2*sig2)
      cc <- 2/sig2 + psi^2/(t2*sig2)
      vv <- tapply(bb2, 1:N, rgig, n=1, lambda=1/2, psi=cc)
      vv[vv<lower] <- lower
    
    ## variables in SAR prior
      #tau2
      bbb <- t(th)%*%iQ%*%th/(2*sig2) + b_tau
      tau2 <- min(max(rinvgamma(1, n/2+a_tau, bbb),10^(-5)), 10^5)   
      
      # rho (random walk MH)
      new.rho <- runif(1,rho-delta,rho+delta)
      if(new.rho < 0){ new.rho <- abs(new.rho) }else if(new.rho > 1){ new.rho <- 2-new.rho }
      new.IW <- as.matrix(diag(1,n)-new.rho*W)
      if(sparse==T){
        new.TIW <- as(t(new.IW), "sparseMatrix")
        new.IW <- as(new.IW, "sparseMatrix")
        new.iQ <- as.matrix(new.TIW%*%(new.IW))
        new.detiQ <- Eigen_det(asSparseMatrix(new.iQ))
      }else{
        new.iQ <- t(new.IW)%*%(new.IW)
        new.detiQ <- det(new.iQ)
      }
      log.den <- -0.5*log(1/new.detiQ)-t(th)%*%new.iQ%*%th/(2*tau2*sig2)  # denominator: log(full)
      log.num <- -0.5*log(1/detiQ)-t(th)%*%iQ%*%th/(2*tau2*sig2) # numerator: log(full)
      prob.r <- min(exp(log.den - log.num), 1)
      if(runif(1)<prob.r){
        rho <- new.rho
        ac.rho <- ac.rho + 1
        iQ <- new.iQ  # update matrix Q
        detiQ <- new.detiQ
      }
      rho.pos[k] <- rho
      
    # print 
    if(TorF==T){
      if(round(k/n.print)==(k/n.print)){  print(k) }
      #if(round(k/n.print)==(k/n.print)){ print(paste0("[",k,"] ac.rho:",ac.rho/k)) }
    }
  }
  #om <- 1:bn
  #th.pos <- th.pos[-om,]
  return(list(th.pos=th.pos, sig.pos=sig.pos, rho.pos=rho.pos, ac.rho=ac.rho/mc))
}







# -------------------------------------------------- #
#         Bayesian spatial quantile smoothing        #
#                   under GP prior                   #
# -------------------------------------------------- #

BSQS.GP <- function(Coordinate, Y, pp=0.2, mc=1000, delta=0.3, coordinatewise=F, n.print=F){
  
  if(n.print==F){TorF <- F}else{TorF <- T}
  if(dim(Coordinate)[1]!=length(Y) | dim(Coordinate)[2]!=2){stop("The dimention of Coordinate is N*2, where N=length(Y)")}
  
  # dataset
  X1 <- Coordinate[,1]
  X2 <- Coordinate[,2]
  Data <- data.frame(cbind(Coordinate,Y))
  colnames(Data) <- c("X1","X2","Y")
  Data <- Data[order(Data$X1),]
  X1 <- Data$X1
  X2 <- Data$X2
  X <- as.data.frame(cbind(X1,X2))
  Y <- Data$Y
  location <- unique(X)
  n <- dim(location)[1]
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- meany <- rep(NA, n)
  num <- rep(0, n+1)
  Xvec <- Yvec <- NULL
  for (i in 1:n) {
    data <- subset(Data, (Data$X1 == location[i,1] & Data$X2 == location[i,2]))
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
    meany[i] <- mean(data$Y)
    Xvec <- c(Xvec, rep(i, Ni[i]))
    Yvec <- c(Yvec, data$Y)
  }
  Y <- Yvec

  ## settings
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  
  # initial values for matrix Q
  rho <- 0.5
  Q1 <- matrix(0, n,n)
  for (i in 1:n){
    for (j in i:n){
      #Q1[j,i] <- Q1[i,j] <- exp(-norm(location[i,]-location[j,],"2")^2)
      Q1[j,i] <- Q1[i,j] <- exp(-norm(location[i,]-location[j,],"2"))
    }
  }
  make_Q <- function(rho){
    Q <- Q1^(1/(2*rho))
    return(Q)
  }
  Q <- make_Q(rho)
  iQ <- solve(Q)
  ac.rho <- 0
  
  ## MCMC box
  th.pos <- matrix(NA, mc, n)
  sig.pos <- c()
  rho.pos <- c()
  
  
  ## initial value
  th <- meany
  sig2 <- 1
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  tau2 <- 1   # 
  a_tau <- b_tau <- 0.1
  vv <- rep(1, N)   # latent parameter in asymmetric Laplace,Exp(sig2)
  IGa <- 1   # paramater of sig2 prior
  IGb <- 1   # paramater of sig2 prior
  

  ## MCMC replications
  for(k in 1:mc){
    # theta
    if(coordinatewise==T){
      TDWD <- iQ/tau2
      A1 <- as.vector(sapply(split(vv^(-1), Xvec), sum))/t2 + diag(TDWD)
      B1 <- as.vector(sapply(split((Y-psi*vv)/vv, Xvec), sum))/t2
      for (i in 1:n) {
        ai <- A1[i]
        bi <- B1[i] - t(TDWD[i,-i])%*%th[-i]
        th[i] <- rnorm(1, bi/ai, sqrt(sig2/ai))
      }
      th.pos[k,] <- th
    }else{
      cc <- as.vector(sapply(split(vv^(-1), Xvec), sum))/t2
      A <- solve(iQ/tau2 + diag(cc) )
      A <- (t(A) + A)/2
      B <- as.vector(sapply(split((Y-psi*vv)/vv, Xvec), sum))/t2
      th <- mvrnorm(1, A%*%B, sig2*A)
      th.pos[k,] <- th
    }
    
    ## asymetric Laplace 
    # sig2
    bb <- rep(0, N)
    for (i in 2:(n+1)) {
      bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-th[i-1]
    }
    cc <- sum(as.vector((bb-psi*vv)^2)/vv)/(2*t2) + sum(vv) + t(th)%*%iQ%*%th/tau2 + IGb
    sig2 <- rinvgamma(1, (n+3*N)/2+IGa, cc)   
    sig.pos[k] <- sig2
    
    # vv
    bb2 <- bb^2/(t2*sig2)
    cc <- 2/sig2 + psi^2/(t2*sig2)
    vv <- tapply(bb2, 1:N, rgig, n=1, lambda=1/2, psi=cc)
    vv[vv<lower] <- lower
    
    ## variables in GP prior
    #tau2
    bbb <- t(th)%*%iQ%*%th/(2*sig2) + b_tau
    tau2 <- min(max(rinvgamma(1, n/2+a_tau, bbb),10^(-5)), 10^5) 
    
    # rho (random walk MH)
    new.rho <- runif(1,rho-delta,rho+delta)
    if(new.rho < 0){ new.rho <- abs(new.rho) }else if(new.rho > 1){ new.rho <- 2-new.rho }
    new.Q <- make_Q(new.rho)
    new.iQ <- solve(new.Q)
    log.den <- -0.5*log(det(new.Q))-t(th)%*%new.iQ%*%th/(2*tau2*sig2)  # denominator: log(full)
    log.num <- -0.5*log(det(Q))-t(th)%*%iQ%*%th/(2*tau2*sig2) # numerator: log(full)
    prob.r <- min(exp(log.den - log.num), 1)
    if(runif(1)<prob.r){
      rho <- new.rho
      ac.rho <- ac.rho + 1
      iQ <- new.iQ  # update matrix Q
      Q <- new.Q
    }
    rho.pos[k] <- rho
    
    # print 
    if(TorF==T){
      if(round(k/n.print)==(k/n.print)){  print(k) }
      #if(round(k/n.print)==(k/n.print)){ print(paste0("[",k,"] ac.rho:",ac.rho/k)) }
    }
  }
  #om <- 1:bn
  #th.pos <- th.pos[-om,]
  return(list(th.pos=th.pos, sig.pos=sig.pos, rho.pos=rho.pos, ac.rho=ac.rho/mc))
}



