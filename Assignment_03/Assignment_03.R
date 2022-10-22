library(tidyverse)

# m
m <- function(x){
  return(sin(3/(x+0.3)))
}

# Second derivative of m
m_prime_prime <- function(x) {
  return(((6*x+1.8)*cos(3/(x + 0.3)) - 9*sin(3/(x+0.3)))/(x+0.3)^4)
}

# Create Location values
locations <- function(n, shape1, shape2){
  return(rbeta(n, shape1, shape2))
}

# Create Response values with standard Gaussian distributed error
resp_val <- function(x, m){
  return(sapply(x, function(z){m(z)+rnorm(1, mean = 0, sd = 1)}))
}

# Function computing h_opt for Gaussian Kernel with Integral over squared kernel 
# density equal to one and integral of second moment with kernel density equal to 1/(2*pi)
h_opt <- function(x, n, m_2){
  return(n^(-1/5)*(1/(2*sqrt(pi)*(m_2(x)*1)^2))^(1/5))
}

# Set the stage for Plots in R Markdown
size <- c(40, 200)
shape1 <- c(1.5, 20.66)
shape2 <- c(100.8, 5)
plot_headings <- list()
for (s in size){
  for (s1 in shape1){
    for (s2 in range(shape2)){
      plot_headings <- append(plot_headings, sprintf("Size = %.2f, Shape1 = %.2f, Shape2 = %.2f", s, s1, s2))
    }
  }
}

# Plots
for (s in size){
  for (s1 in shape1){
    for (s2 in range(shape2)){
      x <- sort(locations(s, s1, s2))
      y <- resp_val(x, m)
      par(mfrow=c(3, 1))
      plot(x, y, xlab = "Locations X", ylab = "Response Values Y")
      lines(x, m(x))
      plot(x, h_opt(x, s, m_prime_prime), "l", xlab = "Locations X", ylab = expression(h[opt](X)), col = "blue")
      plot(x, dbeta(x, shape1 = s1, shape2 = s2), "l", xlab = "Locations X", ylab = expr("pmf"~beta(paste(!!s1, ", " ,!!s2))), col = "red")
    }
  }
}



