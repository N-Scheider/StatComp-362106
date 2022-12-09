# Imported from Assignment
dposterior <- function(x, y, scale = 1.2, sd = .6) {
  # x: evaluation points of the density
  # y: observation Y=y (length 1),
  # scale: scale parameter of Weibull prior (shape=2 fixed)
  # sd: standard deviation of Gaussian error (mean=0 fixed) 
  a <- 1/2*1/sd^2
  c <- 1/scale^2
  erf <- function(x) 2*pnorm(x*sqrt(2)) - 1
  k <- ifelse(x >= 0, x * exp( -a * (x-y)^2 - c*x^2 ), 0) 
  n <- exp(-a*(y^2))*(sqrt(pi)*a*y*exp(a^2*y^2/(a+c))*(erf(a*y/sqrt(a+c))+1)+sqrt(a+c))/(2*(a+c)^(3/2))
  return(k/n)
}

MetropolisHasting <- function(N, y){
  MC <- rep(0, N+1)
  acc <- rep(0, N)
  MC[1] <- y
  for (n in 1:N) {
    X <- rnorm(1, mean=y, sd=0.6)
    acc[n] <- min(1, dweibull(X, shape=2, scale=1.2)/dweibull(MC[n], shape=2, scale=1.2)*dnorm(MC[n], mean=y, sd=0.6)/dnorm(X, mean=y, sd=0.6))
    u <- runif(1, min=0, max=1)
    MC[n+1] <- X*(u<=acc[n]) + MC[n]*(u>acc[n])
  }
  return(list(MC, acc))
}


N <- 1000
burn_in <- 1000

MC_a <- MetropolisHasting(burn_in+N, 0.5)
acc_a <- mean(MC_a[[2]][-c(1:burn_in)])
MC_a <- MC_a[[1]][-c(1:burn_in)]
MC_b <- MetropolisHasting(burn_in+N, 1)
acc_b <- mean(MC_b[[2]][-c(1:burn_in)])
MC_b <- MC_b[[1]][-c(1:burn_in)]
MC_c <- MetropolisHasting(burn_in+N, 2)
acc_c <- mean(MC_c[[2]][-c(1:burn_in)])
MC_c <- MC_c[[1]][-c(1:burn_in)]

#dev.off()
#plot(1:(N+1), MC_a, "l", col=1, xlab="nth step", ylab="state")
#lines(1:(N+1), MC_b, col=2)
#lines(1:(N+1), MC_c, col=3)
#title(main="MC  with Metropolis Hasting")
#legend("topright", legend=c("y=0.5", "y=1", "y=2"), col=c(1,2,3), pch=16)

par(mfrow=c(1, 3))
# given observation y=0.5
hist(MC_a, ylim=c(0,1.3), probability=T, xlab="state, y=0.5", main="", breaks = seq(min(MC_a), max(MC_a), length.out = 25))
lines(seq(min(MC_a), max(MC_a), 0.1), dposterior(seq(min(MC_a), max(MC_a), 0.1), 0.5), col=1)
lines(seq(min(MC_a), max(MC_a), 0.1), dnorm(seq(min(MC_a), max(MC_a), 0.1), 0.5, 0.6), col=2)
# given observation y=1
hist(MC_b, ylim=c(0,1.3), probability=T, xlab="state, y=1", main="", breaks = seq(min(MC_b), max(MC_b), length.out = 25))
lines(seq(min(MC_b), max(MC_b), 0.1), dposterior(seq(min(MC_b), max(MC_b), 0.1), 1), col=1)
lines(seq(min(MC_b), max(MC_b), 0.1), dnorm(seq(min(MC_b), max(MC_b), 0.1), 1, 0.6), col=2)
# given observation y=2
hist(MC_c, ylim=c(0,1.3), probability=T, xlab="state, y=2", main="", breaks = seq(min(MC_c), max(MC_c), length.out = 25))
lines(seq(min(MC_c), max(MC_c), 0.1), dposterior(seq(min(MC_c), max(MC_c), 0.1), 2), col=1)
lines(seq(min(MC_c), max(MC_c), 0.1), dnorm(seq(min(MC_c), max(MC_c), 0.1), 2, 0.6), col=2)
title("Histogram, Posteriori and Proposal densities", line = - 2, outer = TRUE)

print(paste0("Acceptance rate for y=0.5: ", acc_a))
print(paste0("Acceptance rate for y=1: ", acc_b))
print(paste0("Acceptance rate for y=2: ", acc_c))


