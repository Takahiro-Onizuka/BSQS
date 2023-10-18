# --------------------------------------------------------- #
#             noise distribution (flat, four, wave)         #
# --------------------------------------------------------- #


# ------------------------------- #
#             flat noise          #
# ------------------------------- #

flatnoise <- function(th0, kappa=1, pp.set=c(0.1, 0.3, 0.5, 0.7, 0.9) , N=100, seed=123){
  dim1 <-  dim2 <- 10 
  p <- dim1*dim2
  n <- dim1*dim2
  #th0 <- matrix(0,dim1,dim2)
  L <- length(pp.set)
  TT <- th0
  Tquan <- matrix(NA,n,L)
  for (i in 1:L) {
    TT <- th0 + qnorm(pp.set[i],0,kappa)
    Tquan[,i] <- TT
  }
  
  # noise
  set.seed(seed)  # seed
  Y <- matrix(NA, n, N)
  YY <- th0
  for(j in 1:N){
    YY <- th0 + rnorm(n,0,kappa)
    Y[,j] <- as.numeric(YY)
  }
  
  return(list(Tquan=Tquan, Y=Y))
}



# ------------------------------- #
#             four noise          #
# ------------------------------- #

fournoise <- function(th0, kappa=c(1,2,3), pp.set=c(0.1, 0.3, 0.5, 0.7, 0.9) , N=100, seed=123){
  dim1 <-  dim2 <- 10 
  p <- dim1*dim2
  n <- dim1*dim2
  #th0 <- matrix(0,dim1,dim2)
  L <- length(pp.set)
  TT <- th0
  Tquan <- matrix(NA,n,L)
  for (i in 1:L) {
    TT[row(TT)<=5 & col(TT) <=5] <- th0[row(TT)<=5 & col(TT) <=5] + qnorm(pp.set[i],0,kappa[1])
    TT[row(TT)>5 & col(TT) <=5] <- th0[row(TT)>5 & col(TT) <=5] + qnorm(pp.set[i],0,kappa[2])
    TT[row(TT)<=5 & col(TT) >5] <- th0[row(TT)<=5 & col(TT) >5] + qnorm(pp.set[i],0,kappa[2])
    TT[row(TT)>5 & col(TT) >5] <- th0[row(TT)>5 & col(TT) >5] + qnorm(pp.set[i],0,kappa[3])
    Tquan[,i] <- TT
  }
  # noise
  set.seed(seed)  # seed
  Y <- matrix(NA, n, N)
  YY <- th0
  for(j in 1:N){
    YY[row(YY)<=5 & col(YY) <=5] <- th0[row(YY)<=5 & col(YY) <=5] + rnorm(25,0,kappa[1])
    YY[row(YY)>5 & col(YY) <=5] <- th0[row(YY)>5 & col(YY) <=5] + rnorm(25,0,kappa[2])
    YY[row(YY)<=5 & col(YY) >5] <- th0[row(YY)<=5 & col(YY) >5] + rnorm(25,0,kappa[2])
    YY[row(YY)>5 & col(YY) >5] <- th0[row(YY)>5 & col(YY) >5] + rnorm(25,0,kappa[3])
    Y[,j] <- as.numeric(YY)
  }
  
  return(list(Tquan=Tquan, Y=Y))
}




# ------------------------------- #
#             wave noise          #
# ------------------------------- #

wavenoise <- function(th0, kappa=c(1,2,3,4), pp.set=c(0.1, 0.3, 0.5, 0.7, 0.9) , N=100, seed=123){
  dim1 <-  dim2 <- 10 
  p <- dim1*dim2
  n <- dim1*dim2
  #th0 <- matrix(0,dim1,dim2)
  L <- length(pp.set)
  TT <- th0
  Tquan <- matrix(NA,n,L)
  for (i in 1:L) {
    TT[row(TT)<=4 & col(TT) <=4] <- th0[row(TT)<=4 & col(TT) <=4] + qnorm(pp.set[i],0,kappa[1])
    TT[(row(TT)>4 &row(TT)<=6) | (col(TT)>4 &col(TT) <=6)] <- th0[(row(TT)>4 &row(TT)<=6) | (col(TT)>4 &col(TT) <=6)] + qnorm(pp.set[i],0,kappa[2])
    TT[(row(TT)>6 &row(TT)<=8) | (col(TT)>6 &col(TT) <=8)] <- th0[(row(TT)>6 &row(TT)<=8) | (col(TT)>6 &col(TT) <=8)] + qnorm(pp.set[i],0,kappa[3])
    TT[row(TT)>8 | col(TT) >8] <- th0[row(TT)>8 & col(TT) >8] + qnorm(pp.set[i],0,kappa[4])
    Tquan[,i] <- TT
  }
  
  # noise
  set.seed(seed)  # seed
  Y <- matrix(NA, n, N)
  YY <- th0
  for(j in 1:N){
    YY[row(YY)<=4 & col(YY) <=4] <- th0[row(YY)<=4 & col(YY) <=4] + rnorm(16,0,kappa[1])
    YY[(row(YY)>4 &row(YY)<=6 & col(YY)<=6) | (col(YY)>4 &col(YY) <=6& row(YY)<=6)] <- th0[(row(YY)>4 &row(YY)<=6& col(YY)<=6) | (col(YY)>4 &col(YY)<=6& row(YY)<=6)] + rnorm(20,0,kappa[2])
    YY[(row(YY)>6 &row(YY)<=8& col(YY)<= 8) | (col(YY)>6 &col(YY) <=8& row(YY)<=8)] <- th0[(row(YY)>6 &row(YY)<=8& col(YY)<= 8) | (col(YY)>6 &col(YY) <=8& row(YY)<=8)] + rnorm(28,0,kappa[3])
    YY[row(YY)>8 | col(YY) >8] <- th0[row(YY)>8 & col(YY) >8] + rnorm(36,0,kappa[4])
    Y[,j] <- as.numeric(YY)
  }
  
  return(list(Tquan=Tquan, Y=Y))
}



