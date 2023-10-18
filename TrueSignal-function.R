# ------------------------------------------------- #
#     true signals (two block, exponantial, VS)     #
# ------------------------------------------------- #
library(MASS)

# ----------------------- #
#        two block        #
# ----------------------- #


twoBlock <- function(kappa=5){
  dim1 <-  dim2 <- 10 
  p <- dim1*dim2
  n <- dim1*dim2
  th0 <- matrix(0,dim1,dim2)
  th0[(row(th0)-dim1/2)^2 + (col(th0)-dim2/2)^2 <= 4] =  kappa
  th0[(row(th0)-1-dim1/2)^2 + (col(th0)-dim2/2)^2 <= 4] = kappa
  th0[(row(th0)-dim1/2)^2 + (col(th0)-1-dim2/2)^2 <= 4] = kappa
  th0[(row(th0)-1-dim1/2)^2 + (col(th0)-1-dim2/2)^2 <= 4] = kappa
  
  return(th0)
}


# ----------------------- #
#       exponential       #
# ----------------------- #


expBlock <- function(L=5, mu=c(5.5, 5.5), Sigma=matrix(c(3,0,0,3), 2, 2)){
  dim1 <-  dim2 <- 10 
  p <- dim1*dim2
  n <- dim1*dim2
  #L <- 90; mu <- c(5.5, 5.5); Sigma <- matrix(c(2.5,0,0,2.5), 2,2)
  #L*dmvnorm(c(i,j), mu, Sigma)
  th0 <- matrix(0,dim1,dim2)
  exp2 <- function(L, x, mun, Sigma){
    t(x-mu)%*%Sigma%*%(x-mu)
  }
  for(i in 1: dim1){
    for(j in 1:dim2){
      th0[row(th0)==i & col(th0)==j] <- round(L*exp(-t(c(i,j)-mu)%*%solve(Sigma)%*%(c(i,j)-mu)/2), digits = 1)
    }
  }
  return(th0)
}


