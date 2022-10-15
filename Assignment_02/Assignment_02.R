# rmixnorm and dmixnorm imported from the Module 09_GaussianMixture.html

# returns RV for a Gaussian mixture
rmixnorm <- function(N, mu1, mu2, sigma1, sigma2, tau){
  ind <- I(runif(N) > tau)
  X <- rep(0,N)
  X[ind] <- rnorm(sum(ind), mu1, sigma1)
  X[!ind] <- rnorm(sum(!ind), mu2, sigma2)
  return(X)
}

# theoretical pdf for a Gaussian mixture
dmixnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dnorm(x,mu1,sigma1) + tau*dnorm(x,mu2,sigma2)
  return(y)
}

n = 200

bw = seq(0.1, 0.9, 0.05)

gau_err <- matrix(nrow = n, ncol = length(bw))
epane_err <- matrix(nrow = n, ncol = length(bw))
rect_err <- matrix(nrow = n, ncol = length(bw))

for (i in 1:n) {
  t = 0.25
  # Generate 100 samples from the Gaussian mixture
  rv <- rmixnorm(100, 0, 0, 1, 1, t)
  # density estimation for different Kernels with Bandwiths 0.1, 0.15, 0.2
  for (k in 1:length(bw)){
    # error of Gaussian Kernel
    gaussian <- density(rv, kernel="gaussian", bw=bw[k])
    gau_err[i, k] <- norm(sapply(gaussian$x, function(z) dmixnorm(z, 0, 0, 1, 1, t))-gaussian$y, type = c("2"))
    # error of Epanechnikov Kernel
    epane <-density(rv, kernel="epanechnikov", bw=bw[k])
    epane_err[i, k] <- norm(sapply(epane$x, function(z) dmixnorm(z, 0, 0, 1, 1, t))-epane$y, type = c("2"))
    # error of Rectangular Kernel
    rect <- density(rv, kernel="rectangular", bw=bw[k])
    rect_err[i, k] <- norm(sapply(rect$x, function(z) dmixnorm(z, 0, 0, 1, 1, t))-rect$y, type = c("2"))
  }
}

y1 = rep(0, length(bw))
y2 = rep(0, length(bw))
y3 = rep(0, length(bw))

for (k in 1:length(bw)){
  y1[k] = mean(gau_err[,k])
  y2[k] = mean(epane_err[,k])
  y3[k] = mean(rect_err[,k])
}

par(mfrow=c(1,3))
plot(bw, y1, col = "red", type = "l", xlab = "Bandwidth", ylab = "Err Gaussian Kernel")
plot(bw, y2, col = "blue", type = "l", xlab = "Bandwidth", ylab = "Err Epanechnikov Kernel")
plot(bw, y3, col = "green", type = "l", xlab = "Bandwidth", ylab = "Err Rectangular Kernel")
mtext("Error plots of Kernels to Density", side = 3, line = - 2, outer = TRUE)


