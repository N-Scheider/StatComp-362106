# Imported from Module 9
rmixnorm <- function(N, mu1, mu2, sigma1, sigma2, tau){
  ind <- I(runif(N) > tau)
  X <- rep(0,N)
  X[ind] <- rnorm(sum(ind), mu1, sigma1)
  X[!ind] <- rnorm(sum(!ind), mu2, sigma2)
  return(X)
}

dmixnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dnorm(x,mu1,sigma1) + tau*dnorm(x,mu2,sigma2)
  return(y)
}


# Returns expected value of given parameters for gaussian mixed distribution
ExpectationStep <- function(X, mu1, mu2, sigma1, sigma2, tau){
  N <- length(X)
  sum_p <- sum(sapply(X, function(x){p(x, mu1, mu2, sigma1, sigma2, tau)}))
  
  sum1 <- log(1-tau)*(N-sum_p)
  sum2 <- log(tau)*(sum_p)
  sum3 <- sum(sapply(X, function(x){(1-p(x, mu1, mu2, sigma1, sigma2, tau))*log(dnorm((x-mu1)/sigma1))}))
  sum4 <- sum(sapply(X, function(x){p(x, mu1, mu2, sigma1, sigma2, tau)*log(dnorm((x-mu2)/sigma2))}))
  
  return(sum1+sum2+sum3+sum4)
}

# Returns supporting function p for EM of mixed gaussian distribution
p <- function(x, mu1, mu2, sigma1, sigma2, tau){
  return(dnorm((x-mu2)/sigma2)*tau/dmixnorm(x, mu1, mu2, sigma1, sigma2, tau))
}

# Returns optimization of parameters for EM of mixed gaussian distribution
MaximizationStep <- function(X, mu1, mu2, sigma1, sigma2, tau){
  N <- length(X)
  sum_p <- sum(sapply(X, function(x){p(x, mu1, mu2, sigma1, sigma2, tau)}))
  sum_p_min <- sum(sapply(X, function(x){1-p(x, mu1, mu2, sigma1, sigma2, tau)}))
  
  tau_up <- 1/((N-sum_p)/sum_p+1)
  mu1_up <- sum(sapply(X, function(x){(1-p(x, mu1, mu2, sigma1, sigma2, tau))*x}))/sum_p_min
  sigma1_up <- sqrt(2*sqrt(2*pi)*sum(sapply(X, function(x){(1-p(x, mu1, mu2, sigma1, sigma2, tau))*(x-mu1)}))/sum_p_min)
  mu2_up <- sum(sapply(X, function(x){p(x, mu1, mu2, sigma1, sigma2, tau)*x}))/sum_p
  sigma2_up <- sqrt(2*sqrt(2*pi)*sum(sapply(X, function(x){(p(x, mu1, mu2, sigma1, sigma2, tau))*(x-mu2)}))/sum_p)
  
  return(c(mu1_up, mu2_up, sigma1_up, sigma2_up, tau_up))
}



# Test algorithm

N <- 100
X <- rmixnorm(N, 1, 3, 4, 4, 0.5)

#Initialization

mu1 <- 1; mu2 <- 2
sigma1 <- 2; sigma2 <- 5
tau <- 0.3

a <- ExpectationStep(X, mu1, mu2, sigma1, sigma2, tau)


for (i in 1:20) {
  browser()
  b <- MaximizationStep(mu1, mu2, sigma1, sigma2, tau)
  c <- ExpectationStep(X, b[1], b[2], b[3], b[4], b[5])
}

print(b)
print(c)