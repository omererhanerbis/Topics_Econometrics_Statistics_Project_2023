## GMM estimation each step is first provided. Monte Carlo Simulations follow later.

rm(list = ls())
set.seed(1000)


library(gmm)

# Parameters
n = 1000
true_p <- 1 / 2
dependent_mean <- 10

# DGP
X <- rbinom(n, 1, true_p)
U <- runif(n, 0, 1)
t <- as.numeric(true_p >= U)
Y <- rnorm(n, dependent_mean, 1)
Y_obs <- Y * t

var_mat <- as.matrix(data.frame(Y, t, X))

# MM_tw
g_tw <- function(beta, x) {
  m1 <- x[,1] * x[,2] / (1 / 2) - beta
  f <- cbind(m1)
  return(f)
}

# MM_ew
g_ew <- function(beta, x) {
  m1 <- x[,1] * x[,2] / (1 / 2) - beta
  m2 <- x[,3] * (x[,2] - (1 / 2))
  m3 <- (1 - x[,3]) * (x[,2] - (1 / 2))
  f <- cbind(m1, m2, m3)
  return(f)
}


gmm(g_tw, var_mat, c(-1000, 1000), optfct = ("optimize"))$vcov
gmm(g_ew, var_mat, c(-1000, 1000), optfct = ("optimize"))$vcov


gmm_tw <- gmm(g_tw, var_mat, c(-1000, 1000), optfct = ("optimize"))
gmm_tw$coefficients
gmm_tw$vcov
gmm_ew <- gmm(g_ew, var_mat, c(-1000, 1000), optfct = ("optimize"))
gmm_ew$coefficients
gmm_ew$vcov



# Monte Carlo Simulations for GMM

rm(list = ls())
set.seed(200)



gmm_monte_carlo <- function(n_sims = 1000, n = 1000, true_p = 1 / 2, dependent_mean = 10)
{
  tw_estimates <- c()
  ew_estimates <- c()
  for(i in 1:n_sims)
  {
    # DGP
    X <- rbinom(n, 1, true_p)
    U <- runif(n, 0, 1)
    t <- as.numeric(true_p >= U)
    Y <- rnorm(n, dependent_mean, 1)
    Y_obs <- Y * t
    
    var_mat <- as.matrix(data.frame(Y, t, X))
    
    # MM_tw
    g_tw <- function(beta, x) {
      m1 <- x[,1] * x[,2] / (1 / 2) - beta
      f <- cbind(m1)
      return(f)
    }
    # True weights estimator
    sd_tw <- gmm(g_tw, var_mat, c(-1000, 1000), optfct = ("optimize"))$coefficients
    
    
    # MM_ew
    g_ew <- function(beta, x) {
      m1 <- x[,1] * x[,2] / (1 / 2) - beta
      m2 <- x[,3] * (x[,2] - (1 / 2))
      m3 <- (1 - x[,3]) * (x[,2] - (1 / 2))
      f <- cbind(m1, m2, m3)
      return(f)
    }
    # Estimated weights estimator
    sd_ew <- gmm(g_ew, var_mat, c(-1000, 1000), optfct = ("optimize"))$coefficients
    
    
    # Saving estimates to holders
    tw_estimates <- c(tw_estimates, sd_tw)
    ew_estimates <- c(ew_estimates, sd_ew)
    sd_mat <- data.frame(sd_tw_est = tw_estimates, sd_ew_est = ew_estimates)
  }
  return(sd_mat)
}

# Monte Carlo simulations
sds <- gmm_monte_carlo()



#Means
mean(sds[,1])
mean(sds[,2])
# Variances
sd(sds[,1])^2
sd(sds[,2])^2


hist(sds[,1], breaks = 20, main = "True Weights Estimation Distribution", xlab = "beta_0 values")
hist(sds[,2], breaks = 30, main = "True Weights Estimation Distribution", xlab = "beta_0 values")


