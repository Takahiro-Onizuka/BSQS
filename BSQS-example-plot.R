library(igraph)

source("BSQS-function.R")
source("TrueSignal-function.R")
source("NoiseDistribution-function.R")

# MCMC setting
NN <- 5 ; mcs <- 7500; sss <- 10 ; thin <- seq(sss,mcs, by=sss)
pp.set <- 0.5  # quantile level

# 10*10 lattice graph for simulation study
dim1 <-  dim2 <- 10 
glattice <- make_lattice(c(dim1, dim2))
el <- as_edgelist(glattice, names = TRUE)
V <- length(as_ids(V(glattice))) #number of vertex
E <- length(as_ids(E(glattice))) #number of edg
n <- V

# Define the graph difference operator
im <- matrix(0,E,V) 
for (i in 1:E) {
  for (j in 1:V) {
    if(el[i,1]==j){ im[i,j] <- 1 }
    if(el[i,2]==j){ im[i,j] <- -1 }
  }
}
D <- t(im)%*%im

# Define the contingency matrix
W <- as_adjacency_matrix(glattice)
sW <- matrix(0,V,V)
for(i in 1:n){ sW[i,] <- W[i,]/sum(W[i,]) }

# Define the observation point
xx2 <- rep(c(1:10),10)
xx1 <- sort(xx2)
coords <- cbind(xx1,xx2)
for (i in 1:(NN-1)) {
  coords <- rbind(coords, cbind(xx1,xx2))
}


# Data-generating
x1 <- 1:dim1
x2 <- 1:dim2
th0 <- expBlock(L=5, mu=c(5.5, 5.5), Sigma=matrix(c(3,0,0,3), 2, 2))
TquanY <- flatnoise(th0, kappa=1, pp.set=pp.set, N=5, seed=123)
y <- as.vector(TquanY$Y)
x <- rep(1:100,NN)


# Fitting of TF-HS model
set.seed(1)
sample <- BSQTF.HS(x, y, D, pp=pp.set, mc=mcs, n.print=1000)
sm <- apply(sample[thin, ], 2, mean) 


# Fitting of SAR model
set.seed(1)
sample <- BSQS.SAR(X=x, Y=y, W=sW, pp=pp.set, mc=mcs, delta=0.2, n.print=1000)
smSAR <- apply(sample$th.pos[thin, ], 2, mean) 


# Fitting of GP model
set.seed(1)
sample <- BSQS.GP(Coordinate=coords, Y=y, pp=pp.set, mc=mcs, delta=0.2, n.print=1000)
smGP <- apply(sample$th.pos[thin, ], 2, mean) 


## size = 20*4
Zlim <- c(-1,6.5)
x1 <- 1:dim1
x2 <- 1:dim2
par(mfrow=c(1,3))
persp(x1, x2, t(matrix(sm,dim1,dim2)), theta = 30,phi = 15,expand = 0.8, shade = 0.5, border = NA, ticktype = "detailed", 
      zlab="", col = "gray70",ltheta = -120, nticks=6, main="BQTF-HS", cex.lab  = 1, cex.axis = 1, cex.main = 1.5, zlim = Zlim)
persp(x1, x2, t(matrix(smSAR,dim1,dim2)), theta = 30,phi = 15,expand = 0.8, shade = 0.5, border = NA, ticktype = "detailed", 
      zlab="", col = "gray70",ltheta = -120, nticks=6, main="SAR", cex.lab  = 1, cex.axis = 1, cex.main = 1.5, zlim = Zlim)
persp(x1, x2, t(matrix(smGP,dim1,dim2)), theta = 30,phi = 15,expand = 0.8, shade = 0.5, border = NA, ticktype = "detailed", 
      zlab="", col = "gray70",ltheta = -120, nticks=6, main="GP", cex.lab  = 1, cex.axis = 1, cex.main = 1.5, zlim = Zlim)



