count <- rep(0, 4)

for (i in 1:1000) {
  N <- 100
  B <- 1000
  alpha <- 0.05
  
  # Generate RV
  e <- rexp(N, rate = 2)
  mean_e <- mean(e)
  
  # Bootstrap (imported from lecture notes)
  boot_strap <- array(sample(e, N*B, replace=TRUE), c(B, N))
  boot_stat_stud <- sapply(1:B, function(b){sqrt(N)*(mean(boot_strap[b,])-mean_e)/sd(boot_strap[b,])})
  boot_stat_non_stud <- sapply(1:B, function(b){sqrt(N)*(mean(boot_strap[b,])-mean_e)})
  boot_stat_sample <- sapply(1:B, function(b){sqrt(N)*(mean(boot_strap[b,])-mean_e)/sd(e)})
  
  # Calculate CI
  asymp_ci <- mean_e + sd(e)*qnorm(1-alpha)/sqrt(N)
  stud_ci <- mean_e + sd_e*quantile(boot_stat_stud, 1-alpha)/sqrt(N)
  non_stud_ci <- mean_e + quantile(boot_stat_non_stud, 1-alpha)/sqrt(N)
  sample_ci <- mean_e + quantile(boot_stat_sample, 1-alpha)*sd(e)/sqrt(N)
  
  lst <- list(asymp_ci, stud_ci, non_stud_ci, sample_ci)
  
  for (c in 1:4) {
    count[c] <- count[c]+(1/2 <= lst[c])
  }
}

cov <- count/1000

df_cov <- data.frame("Coverage_Ratio"=cov)
rownames(df_cov) <- c("asymptotic", "studentized", "non-studentized", "sample-truth-scaled")

df_cov

