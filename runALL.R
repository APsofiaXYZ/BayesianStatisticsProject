#############

# This script fits spatial point data to the model presented in our report. 
# The data: coordinates for points X, Y with marks Z corresponding to stages of cancer.

# Small sample data could be found in the related repository.
# Please install Rcpp and RcppArmadillo packages. AP.
# The corresponding Cpp file with DMH algorithm and the functions attached.

############

library(RcppArmadillo) 
library(Rcpp)

Rcpp::sourceCpp('DMH.cpp')

s<- read.table('11446_6_Info-201p.txt', header = T)
X<-s$x
Y<-s$y
Z<-as.integer(s$z)

Qnumber<-3
C<-0.1
N<-length(X)

Tau<-sqrt(0.1) # sqrt
M<-1*length(X) # should be multiple of length

# iterations for DMC
Iter <- 50000

# define set of possible types for cells
Marks<-1:Qnumber


### PRIOR DISTRIBUTIONS ###
# set the hyperparameters 

Mu_omega <- 1 # as suggested by Li
Sigma_omega <- 1 # as suggested by Li
Mu_theta <- 0 # as suggested by Li
Sigma_theta <- 1 # as suggested by Li
A_lambda <- 0.001 # as suggested by Li, Gelman
B_lambda <- 0.001
Lambda<-rgamma(1, 1, 1)

# ==== Initialization for chain ====
Omega_start <-rep(1,Qnumber)
Omega_start [1]<-rnorm(1, Mu_omega, 1)

# Theta matrix, symmetric, let's generate from priors 
Theta_start <- matrix(0,Qnumber,Qnumber)
for (l in 1:Qnumber) {
  for (k in l:Qnumber) {
    Theta_start[l,k] <-rnorm(1,Mu_theta,Sigma_theta)
    Theta_start[k,l] <-Theta_start[l,k]
    
  }
}
Theta_start[Qnumber,Qnumber]<-1

# ========

Omega <- Omega_start
Theta <- Theta_start

# d-matrix
Dmatrix <- matrix(0,N,N)
for (l in 1:N) {
  for (k in l:N) {
    Dmatrix[l,k] <-abs(sqrt((X[k]-X[l])^2+(Y[k]-Y[l])^2))
    Dmatrix[k,l] <-Dmatrix[l,k]
  }}


# ======== DMH

final<-DMH(X, Y, Z, Marks, Dmatrix, Iter, M, N, Qnumber, C, Lambda, Mu_omega, Sigma_omega, Tau, Mu_theta, Sigma_theta, A_lambda, B_lambda, Omega_start, Theta_start)

omega_matrix<-final$omega_matrix
acceptomega<-final$accept_omega
accepttheta<-final$accept_theta
theta_matrix<-final$theta_matrix
la_vec<-final$lambda_vec
accept_la<-final$accept_la


save(omega_matrix,acceptomega, accepttheta, theta_matrix, 
     la_vec, accept_la, Tau, Iter, M, N, 
     file = "pat_p100012_201p.RData")
