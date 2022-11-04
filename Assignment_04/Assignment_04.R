library(MASS)
library(locpol)
data(mcycle)
mcycle <- mcycle[mcycle$times <= 40, ]
plot(mcycle$times,mcycle$accel)
p_degree <- c(1, 2, 3)
hgrid <- 3:15

# data frame to capture different errors
df_err <- data.frame(matrix(nrow=length(hgrid), ncol=3))
colnames(df_err) <- c("Degree 1", "Degree 2", "Degree 3")
rownames(df_err) <- c(hgrid)

# length of data set
n <- length(mcycle$times)

# LOO and MSE
for (p in p_degree) {
  for (h in hgrid) {
    mse <- 0
    for (i in 1:n) {
      model <- locpol(formula = accel[-i]~times[-i], data = mcycle, bw = h, deg = p, xeval = mcycle$times[i])
      mse <- mse + (mcycle$accel[i]-model$lpFit[2])^2/n
    }
    df_err[which(rownames(df_err) == h), p] <- mse
  }
}

plot(hgrid, df_err$`Degree 1`, type = 'l', col = 1, ylab = "MSE")
lines(hgrid, df_err$`Degree 2`, col = 2)
lines(hgrid, df_err$`Degree 3`, col = 3)
legend('topleft', legend = c(colnames(df_err)), col = 1:3, pch = 16)

min_min <- which(df_err == min(df_err, na.rm = TRUE), arr.ind = TRUE)
best_model <- locpol(formula = accel~times, data = mcycle, bw = min_min[1], deg = min_min[2], xeval = mcycle$times)
plot(mcycle$times, mcycle$accel)
lines(mcycle$times, best_model$lpFit$accel)
legend('topleft', legend = sprintf("h = %.f, p = %.f, MSE = %.2f", min_min[1], min_min[2], df_err[min_min]))
