# Exercise 1

y <- rnorm(100)
x <- matrix(rnorm(5000), nrow = 100)
linear_model <- lm(y ~ 0 + x)
cat("R squared error in initial model:", summary(linear_model)$r.squared)

x_sign <- x[,which(summary(linear_model)$coefficients[,"Pr(>|t|)"]<0.05)]
linear_model_sign <- lm(y ~ 0 + x_sign)
cat("R squared error in model with sign coef:", summary(linear_model_sign)$r.squared)


# Exercise 2

# Idea 1: Get the n most significant p-values and compute the R squared error of their RV
n = 10
p_limit <- sort(summary(linear_model)$coefficients[,4])[n]
x_sign2 <- x[,which(summary(linear_model)$coefficients[,"Pr(>|t|)"]<p_limit)]
linear_model_sign2 <- lm(y ~ 0 + x_sign2)
cat("R squared error in model with sign coef:", summary(linear_model_sign2)$r.squared)

# Idea 2 (lecture): Append RV
peeking <- function(X_rv_n = 5000, X_rv_rows = 100, ext = 45){
  X_rv <- matrix(rnorm(X_rv_n), nrow = X_rv_rows)
  #  creating function to get T sizes of columns without and with appending new RV
  Tsize_func <- function(vec){
    return(mean(vec)/sd(vec)*sqrt(length(vec)))
  }
  Tsize_func_app <- function(vec){
    vec <- append(vec, rnorm(ext))
    return(mean(vec)/sd(vec)*sqrt(length(vec)))
  }
  # Preserve significant columns and peek insignificant ones
  Xvec_Tstats <- mapply(Tsize_func, asplit(X_rv, 2))
  ind <- which(Xvec_Tstats <= qt(p = 0.975, df = X_rv_rows-1))
  Xvec_Tstats2 <- mapply(Tsize_func_app, asplit(X_rv[,ind], 2))
  Xvec_Tstats <- replace(Xvec_Tstats, ind, Xvec_Tstats2)
  return(mean(I(abs(Xvec_Tstats) > qnorm(0.975))))
}

mean(replicate(100, peeking()))

