library(VGAM)

N <-  1e5
sim <- 50
t_stand <- rep(0, sim)
t_inv <- rep(0, sim)
t_box <- rep(0, sim)
t_rej <- rep(0, sim)
KS_stats <- rep(0, sim)
KS_stats_inv <- rep(0, sim)
KS_stats_box <- rep(0, sim)
KS_stats_rej  <- rep(0, sim)


for (i in 1:sim) {
  # Kolmogorov-Smirnov Statistics
  # s0 <- Sys.time()
  # X <- rnorm(N)
  # F_hat <- ecdf(X)
  # F_hat_eval <- F_hat(sort(X))
  # F_eval <- pnorm(sort(X))
  # KS_stats[i] <- max(abs(F_hat_eval - F_eval))
  # t_stand[i] <- Sys.time() - s0
  
  U <- runif(N)
  U2 <- runif(N)
  
  # Gaussian Samples via inverse transform
  s1 <- Sys.time()
  X_inv <- qnorm(U)
  F_inv_ecdf <- ecdf(X_inv)
  F_inv_eval <- F_inv_ecdf(sort(X_inv))
  KS_stats_inv[i] <- max(abs(F_inv_eval - F_eval))
  t_inv[i] <- Sys.time() - s1
  
  
  
  # Gaussian Samples via Box-Muller transform
  s2 <- Sys.time()
  X_box <- sqrt(-2*log(U))*cos(2*pi*U2)
  F_box_ecdf <- ecdf(X_box)
  F_box_eval <- F_inv_ecdf(sort(X_box))
  KS_stats_box[i] <- max(abs(F_box_eval - F_eval))
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
  F_rej_eval <- F_rej_ecdf(sort(X_rej))
  KS_stats_rej[i] <- max(abs(F_rej_eval - F_eval))
  t_rej[i] <- Sys.time() - s3
}

leg <- c("Inverse", "Box", "Rej")

plot(1:sim, t_inv, type = "l", col = 1)
lines(1:sim, t_box, col = 2)
lines(1:sim, t_rej, col = 3)
legend("topright", legend = leg, col = c(1:3), pch = 15)

plot(1:sim, KS_stats_inv, type = "l", col = 1)
lines(1:sim, KS_stats_box, col = 2)
lines(1:sim, KS_stats_rej, col = 3)
legend("topright", legend = leg, col = c(1:3), pch = 15)


