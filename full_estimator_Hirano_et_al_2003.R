
# This is one step of Monte Carlo simulations to be performed that implements Hirano et al (2003) estimator. Even individual steps are close to improbable to run with high K values

rm(list = ls())
set.seed(321)

# install.packages("devtools")
# library(devtools)
# install_github("ccfang2/dimada")

library(MASS) # for mvnorm
library(polynom)
library(orthopolynom)
library(pracma)
library(dimada)
library(lpSolve)
library(nlme)



# Simulate Data

# Set sample size
n <- 20000

# Set number of covariates
r <- 2

# Set number of terms of basis functions
K <- 1

# Treatment effect size
tau <- 10

# Lower and upper bounds for compact sampling dimensions
lb <- c(0, 0)
ub <- c(pi/4, pi/8)

# DGP
# X defined on cartesian product of compact sets
X_gen <- function(lo_bo=lb, up_bo=ub, size=n, covs=r)
{
  X <- matrix(rep(NA, n*r), ncol = r)
  
  for(i in 1:r)
  {
    X[,i] <- replicate(n, runif(1, lb[i], ub[i]))
  }
  return(X)
}
X <- X_gen()

# rename columns
colnamevector <- c()
for(i in 1:r)
{
  colnamevector <- c(colnamevector, paste("X", i, sep = ""))
}
colnames(X) <- colnamevector
X




# Propensity score function
p_x <- function(X_i){
  f = sin(sum(X_i * c(1, 2)))
  return(f)
}

# True assignment probabilities
true_p <- matrix(rep(NA, n), ncol = 1)
for(i in 1:n)
{
  true_p[i,1] <- p_x(X[i,])
}

# Treatment assignment vector
U <- matrix(runif(n, 0, 1), ncol = 1, byrow = T)
t <- matrix(as.numeric(true_p >= U), ncol = 1, byrow = T)
colnames(t) <- c("T")
t

# Dependent variable in the form of Y = tau * t
Y <- tau * t
colnames(Y) <- c("Y")
Y




# Polynomial basis



power_mat <- poly.gen(X, test.data=NULL, n.basis=K, max.interaction=r, legendre=FALSE)$train


# Initiation values of pi for maximum likelihood estimation
initial_values <- rep(0, ncol(power_mat))




# Negative Likelihood Function Definiton for minimization

negative_likelihood <- function(pi_vec, power_matrix = power_mat, W = t){
  summand <- 0
  for(i in 1:nrow(power_matrix))
  {
    summand <- summand + W[i] * log(exp(sum(power_matrix[i,] * pi_vec))/(1+exp(sum(power_matrix[i,] * pi_vec)))) + (1 - W[i]) * log(exp(sum(1 - power_matrix[i,] * pi_vec))/(1+exp(sum(1 - power_matrix[i,] * pi_vec))))
    # summand <- 
  }
  return(-summand)
}

# Expecting a positive number for negative likelihood minimum value since any value in log is negative for odds ratio domain (0,1)
# negative_likelihood(initial_values)

# SLE
# configure stepmax tuning parameter to enable sufficient step size to estimate maximum to reduce load
SLE <- nlm(negative_likelihood, initial_values)

# estimated nuisance parameter
pi_est <- SLE$estimate
pi_est


# Propensity score estimates
p_est <- exp(as.matrix(power_mat) %*% pi_est ) / (1 + exp(as.matrix(power_mat) %*% pi_est ))
p_est

# True Propensity scores
true_p


# Treatment effects function
tau <- function(p)
{
  summand <- 0
  for(i in 1:n)
  {
    summand <- summand + Y[i,1] * t[i,1] / (p[i,1]) - ( Y[i,1] * (1-t[i,1]) / (1 - p[i,1]) )
  }
  tau <- 1/n * summand
  #ifelse(p == p_est, tau <- 4/3 * tau, tau <- tau)
  return(tau)
}

# Treatment effects estimate using nonparametrically estimated propensity scores by Hirano et al
tau(p_est)
# Treatment effects estimate using true propensity scores
tau(true_p)
