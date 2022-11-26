library(VGAM)

N <-  1e4
sim <- 50

t_inv <- rep(0, sim)
t_box <- rep(0, sim)
t_rej <- rep(0, sim)
KS_stats_inv <- rep(0, sim)
KS_stats_box <- rep(0, sim)
KS_stats_rej  <- rep(0, sim)

for (i in 1:sim) {

  # Gaussian Samples via inverse transform
  s0 <- Sys.time()
  X <- rnorm(N)
  F_inv_ecdf <- ecdf(X)
  F_hat_inv_eval <- F_inv_ecdf(sort(X))
  F_inv_eval <- pnorm(sort(X))
  KS_stats_inv[i] <- max(abs(F_hat_inv_eval - F_inv_eval))
  t_inv[i] <- Sys.time() - s0
  
  # Gaussian Samples via Box-Muller transform
  s2 <- Sys.time()
  U1 <- runif(N)
  U2 <- runif(N)
  X_box <- sqrt(-2*log(U1))*cos(2*pi*U2)
  F_box_ecdf <- ecdf(X_box)
  F_hat_box_eval <- F_box_ecdf(sort(X_box))
  F_box_eval <- pnorm(sort(X_box))
  KS_stats_box[i] <- max(abs(F_hat_box_eval - F_box_eval))
  t_box[i] <- Sys.time() - s2
  
  # Gaussian Samples via rejection sampling
  # c <- sqrt(2/pi)*exp(1/2)
  # xgrid <- seq(-3, 3, 0.01)
  # plot(xgrid, dnorm(xgrid), ylim = c(0, 0.7), type = 'l')
  # lines(xgrid, c*dlaplace(xgrid), col = "blue")
  # legend("topleft", legend = c("proposal", "target"), col = c("green", "black"), pch = 16)
  
  s3 <- Sys.time()
  X_rej <- rep(0,N)
  c <- sqrt(2/pi)*exp(1/2)
  count <- 0
  while(count != N){
    Y <- rlaplace(1)
    Ue <- runif(1)
    if (Ue<1/c*dnorm(Y)/dlaplace(Y)){
      X_rej[count+1] <- Y
      count <- count + 1
    }
  }
  F_rej_ecdf <- ecdf(X_rej)
  F_hat_rej_eval <- F_rej_ecdf(sort(X_rej))
  F_rej_eval <- pnorm(sort(X_rej))
  KS_stats_rej[i] <- max(abs(F_hat_rej_eval - F_rej_eval))
  t_rej[i] <- Sys.time() - s3
}

leg <- c("Inverse", "Box", "Rej")

par(mfrow=c(1, 2))

plot(1:sim, t_inv, ylim = c(0, max(t_rej)), xlab = "simulation", ylab = "time", type = "l", col = 1)
lines(1:sim, t_box, col = 2)
lines(1:sim, t_rej, col = 3)
legend("topright", legend = leg, col = c(1:3), pch = 15)

plot(1:sim, t_inv, type = "l", col = 1, xlab = "simulation", ylab = "time")
lines(1:sim, t_box, col = 2)
legend("topright", legend = leg, col = c(1:2), pch = 15)

dev.off()

plot(1:sim, KS_stats_inv, ylim = c(0, max(KS_stats_box)), , xlab = "simulation", ylab = "Kolmogorov-Smirnov statistic", type = "l", col = 1)
lines(1:sim, KS_stats_box, col = 2)
lines(1:sim, KS_stats_rej, col = 3)
legend("topright", legend = leg, col = c(1:3), pch = 15)


c <- sqrt(2/pi)*exp(1/2)
xgrid <- seq(-3, 3, 0.01)
plot(xgrid, dnorm(xgrid), ylim = c(0, 0.7), type = 'l')
lines(xgrid, c*dlaplace(xgrid), col = "green")
legend("topleft", legend = c("proposal", "target"), col = c("green", "black"), pch = 16)

