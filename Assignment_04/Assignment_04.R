library(MASS)
data(mcycle)
mcycle <- mcycle[mcycle$times <= 40, ]
plot(mcycle$times,mcycle$accel)
