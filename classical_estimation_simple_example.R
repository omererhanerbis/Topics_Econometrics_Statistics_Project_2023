
## If desired, the Monte Carlo Simulations Code is under Step Codes

# Each step Monte Carlo simulations perform

set.seed(1000)


# Simulate data

# Set sample size
n <- 1000

# DGP
X <- rbinom(1000, 1, 0.5)

true_p <- 1 / 2

U <- runif(n, 0, 1)
t <- as.numeric(true_p >= U)

Y <- rnorm(n, 10, 1)
Y_obs <- Y * t

mean(Y_obs, na.rm = T)
mean(Y)

beta_tw <- mean(Y_obs / (1 / 2))

# Population partitions of t * x
N_00 <- sum((X == 0) * (t == 0))
N_10 <- sum((X == 0) * (t == 1))
N_01 <- sum((X == 1) * (t == 0))
N_11 <- sum((X == 1) * (t == 1))

estimated_p <- as.numeric(X == 0) * (N_10 / (N_00 + N_10)) + as.numeric(X == 1) * (N_11 / (N_01 + N_11)) 
beta_ew <- mean(Y_obs / estimated_p)
sd(Y_obs / estimated_p)
sd(Y_obs / true_p)



# Monte Carlo Simulations



# Beta_0s Distributions

rm(list = ls())
set.seed(1000)



monte_carlo <- function(n_sims = 10000, n = 1000, true_p = 1 / 2, dependent_mean = 10)
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
    
    # True weights estimator
    beta_tw <- mean(Y_obs / (1 / 2))
    
    # Selection probabilities 
    N_00 <- sum((X == 0) * (t == 0))
    N_10 <- sum((X == 0) * (t == 1))
    N_01 <- sum((X == 1) * (t == 0))
    N_11 <- sum((X == 1) * (t == 1))
    
    # Estimated weights estimator
    estimated_p <- as.numeric(X == 0) * (N_10 / (N_00 + N_10)) + as.numeric(X == 1) * (N_11 / (N_01 + N_11)) 
    beta_ew <- mean(Y_obs / estimated_p)
    
    # Saving estimates to holders
    tw_estimates <- c(tw_estimates, beta_tw)
    ew_estimates <- c(ew_estimates, beta_ew)
    beta_mat <- data.frame(beta_tw_est = tw_estimates, beta_ew_est = ew_estimates)
  }
  return(beta_mat)
}

# Monte Carlo simulations
betas <- monte_carlo()

#Means
mean(betas[,1])
mean(betas[,2])
# Variances
sd(betas[,1])^2
sd(betas[,2])^2


# Graphs
hist(betas[,1], breaks = 20, main = "True Weights Estimation Distribution", xlab = "beta_0 values")
hist(betas[,2], breaks = 30, main = "Estimated Weights Estimation Distribution", xlab = "beta_0 values")
