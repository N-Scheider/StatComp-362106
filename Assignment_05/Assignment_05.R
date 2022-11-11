# Imported from Module 9
rmixnorm <- function(N, mu1, mu2, sigma1, sigma2, tau){
  ind <- I(runif(N) > tau)
  X <- rep(0,N)
  X[ind] <- rnorm(sum(ind), mu1, sigma1)
  X[!ind] <- rnorm(sum(!ind), mu2, sigma2)
  return(X)
}

dmixnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dnorm(x, mean = mu1, sd = sigma1) + tau*dnorm(x, mean = mu2, sd = sigma2)
  return(y)
}


# Returns expected value of given parameters for gaussian mixed distribution
ExpectationStep <- function(X, mu1, mu2, sigma1, sigma2, tau){
  N <- length(X)
  p <- dnorm(X, mean = mu2, sd = sigma2)*tau/dmixnorm(X, mu1, mu2, sigma1, sigma2, tau)

  sum1 <- log(1-tau)*(N-sum(p))
  sum2 <- log(tau)*sum(p)
  sum3 <- sum((1-p) * log(dnorm(X, mean = mu1, sd = sigma1)))
  sum4 <- sum(p * log(dnorm(X, mean = mu2, sd = sigma2)))
  
  return(sum1+sum2+sum3+sum4)
}

# Returns optimization of parameters for EM of mixed gaussian distribution
MaximizationStep <- function(X, mu1, mu2, sigma1, sigma2, tau){
  N <- length(X)
  p <- dnorm(X, mean = mu2, sd = sigma2)*tau/dmixnorm(X, mu1, mu2, sigma1, sigma2, tau)

  tau_up <- sum(p)/N
  mu1_up <- sum((1-p) * X)/sum(1-p)
  sigma1_up <- sqrt(sum((1-p) * (X-mu1_up)^2)/sum(1-p))
  mu2_up <- sum(p * X)/sum(p)
  sigma2_up <- sqrt(sum(p * (X-mu2_up)^2)/sum(p))
  
  return(c(mu1_up, mu2_up, sigma1_up, sigma2_up, tau_up))
}


# Test algorithm

N <- 1000
X <- rmixnorm(N, 0, 4, 0.5, 2, 0.8)

#Initialization

mu1 <- -1; mu2 <- 6
sigma1 <- 1; sigma2 <- 4
tau <- 0.6

a <- ExpectationStep(X, mu1, mu2, sigma1, sigma2, tau)
b <- MaximizationStep(X, mu1, mu2, sigma1, sigma2, tau)

for (i in 1:1000){
  b <- MaximizationStep(X, b[1], b[2], b[3], b[4], b[5])
  a <- ExpectationStep(X, b[1], b[2], b[3], b[4], b[5])
}

b

xgrid <- seq(from = -4, to = 8, 0.1)
plot(xgrid, dmixnorm(xgrid, 0, 4, 0.5, 2, 0.8), 'l')
lines(xgrid, dmixnorm(xgrid, b[1], b[2], b[3], b[4], b[5]), col = 2)



# Initialization 2

# Initialization

N <- 1000
X <- rmixnorm(N, 0, 4, 0.5, 2, 0.8)

mu1 <- 3; mu2 <- 1
sigma1 <- 3; sigma2 <- 1
tau <- 0.4

a <- ExpectationStep(X, mu1, mu2, sigma1, sigma2, tau)
b <- MaximizationStep(X, mu1, mu2, sigma1, sigma2, tau)

for (i in 1:1000){
  b <- MaximizationStep(X, b[1], b[2], b[3], b[4], b[5])
  a <- ExpectationStep(X, b[1], b[2], b[3], b[4], b[5])
}


print(sprintf("The EM Algorithm approximates the inital parameters as mu1 = %.2f, sigma1 = %.2f, mu2 = %.2f, sigma2 = %.2f, tau = %.2f", b[1], b[2], b[3], b[4], b[5]))

x <- seq(from = -4, to = 8, 0.1)
plot(x, dmixnorm(x, 0, 4, 0.5, 2, 0.8), 'l', xlab = "x", ylab = "density")
lines(x, dmixnorm(x, b[1], b[2], b[3], b[4], b[5]), col = 3)
legend("topright", legend=c("Theoretical pmf", "Approximated pmf"), col = c(1, 2), pch = 16)

