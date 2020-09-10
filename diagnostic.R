library(RcppArmadillo) 
library(Rcpp)
library(spatstat) 
library(coda)


omega1<-omega_matrix
theta1<-theta_matrix
lambda1<-la_vec
acceptomega1<-acceptomega
accepttheta1<-accepttheta
accept_la1<-accept_la

omega2<-omega_matrix
theta2<-theta_matrix
lambda2<-la_vec
acceptomega2<-acceptomega
accepttheta2<-accepttheta
accept_la2<-accept_la

omega3<-omega_matrix
theta3<-theta_matrix
lambda3<-la_vec
acceptomega3<-acceptomega
accepttheta3<-accepttheta
accept_la3<-accept_la


omega1<-omega1[25000:50000,]
theta1<-theta1[25000:50000,]
lambda1<-lambda1[25000:50000]

omega2<-omega2[25000:50000,]
theta2<-theta2[25000:50000,]
lambda2<-lambda2[25000:50000]

omega3<-omega3[25000:50000,]
theta3<-theta3[25000:50000,]
lambda3<-lambda3[25000:50000]

thin = 10
thin_indx = seq(from = 10, 25000, by = 10)

# traceplots
x11()
par(mfrow=c(4,2))
plot(omega1[,1], type = 'l', ylab="Omega 1")
points(omega2[,1], type = 'l', col='red')
points(omega3[,1], type = 'l', col='green')

plot(density(omega1[,1]), type='l', title(main = NULL))
points(density(omega2[,1]), type='l', col='red')
points(density(omega3[,1]), type='l', col='green')

plot(omega1[,2], type = 'l', ylab="Omega 2")
points(omega2[,2], type = 'l', col='red')
points(omega3[,2], type = 'l', col='green')

plot(density(omega1[,2]), type='l', title(main = NULL))
points(density(omega2[,2]), type='l', col='red')
points(density(omega3[,2]), type='l', col='green')

plot(theta1[,1], type = 'l', ylab="Theta 11")
points(theta2[,1], type = 'l', col='red')
points(theta3[,1], type = 'l', col='green')

plot(density(theta1[,1]), type='l')
points(density(theta2[,1]), type='l', col='red')
points(density(theta3[,1]), type='l', col='green')

plot(theta1[,2], type = 'l', ylab="Theta 12")
points(theta2[,2], type = 'l', col='red')
points(theta3[,2], type = 'l', col='green')

plot(density(theta1[,2]), type='l')
points(density(theta2[,2]), type='l', col='red')
points(density(theta3[,2]), type='l', col='green')


plot(theta1[,3], type = 'l', ylab="Theta 13")
points(theta2[,3], type = 'l', col='red')
points(theta3[,3], type = 'l', col='green')

plot(density(theta1[,3]), type='l')
points(density(theta2[,3]), type='l', col='red')
points(density(theta3[,3]), type='l', col='green')

plot(theta1[,4], type = 'l', ylab="Theta 22")
points(theta2[,4], type = 'l', col='red')
points(theta3[,4], type = 'l', col='green')

plot(density(theta1[,4]), type='l')
points(density(theta2[,4]), type='l', col='red')
points(density(theta3[,4]), type='l', col='green')

plot(theta1[,5], type = 'l', ylab="Theta 23")
points(theta2[,5], type = 'l', col='red')
points(theta3[,5], type = 'l', col='green')

plot(density(theta1[,5]), type='l')
points(density(theta2[,5]), type='l', col='red')
points(density(theta3[,5]), type='l', col='green')

plot(lambda1, type = 'l', ylab="Lambda")
points(lambda2, type = 'l', col='red')
points(lambda3, type = 'l', col='green')

plot(density(lambda1), type='l')
points(density(lambda2), type='l', col='red')
points(density(lambda3), type='l', col='green')


# chains
chain1<-unname(cbind(omega1[,1:2],theta1[,1:5], lambda1))
chain2<-unname(cbind(omega2[,1:2],theta2[,1:5], lambda2))
chain3<-unname(cbind(omega3[,1:2],theta3[,1:5], lambda3))

m.list.b <- mcmc.list(as.mcmc(chain1), 
                      as.mcmc(chain2),
                      as.mcmc(chain3))

chain1t <- chain1[thin_indx, ]
chain2t <- chain2[thin_indx, ]
chain3t <- chain3[thin_indx, ]

m.list.bi <- mcmc.list(as.mcmc(chain1t),
                       as.mcmc(chain2t), 
                       as.mcmc(chain3t))

# Gelman-Rubin-Brooks
gelman.diag(m.list.bi)
gelman.plot(m.list.bi)


library(BayesianTools)
correlationPlot(data.frame(chain1t))


# autocorrellation
x11()
par(mfrow = c(3, 3))
acf(omega1[thin_indx, 1], lag.max = 30, lwd = 2)
acf(omega1[thin_indx, 2], lag.max = 30, lwd = 2)
acf(theta1[thin_indx, 1], lag.max = 30, lwd = 2)
acf(theta1[thin_indx, 2], lag.max = 30, lwd = 2)
acf(theta1[thin_indx, 3], lag.max = 30, lwd = 2)
acf(theta1[thin_indx, 4], lag.max = 30, lwd = 2)
acf(theta1[thin_indx, 5], lag.max = 30, lwd = 2)
acf(lambda1[thin_indx], lag.max = 30, lwd = 2)

# ESS
effectiveSize(m.list.bi)

# Heidel
heidel.diag(as.mcmc(chain1t), eps = 0.1, pvalue = 0.05)
heidel.diag(as.mcmc(chain2t), eps = 0.1, pvalue = 0.05)
heidel.diag(as.mcmc(chain3t), eps = 0.1, pvalue = 0.05)

# Geweke
geweke.diag(m.list.bi)

summary(m.list.bi)


